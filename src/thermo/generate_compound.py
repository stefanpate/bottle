# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Institute for Molecular Systems Biology, ETH Zurich.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

'''
KMS CUSTOMIZATION OF EQUILIBRATOR
'''


"""Create Compound objects outside of the compound-cache."""
import logging
from collections import defaultdict
from typing import List, Union

import pandas as pd
from equilibrator_api import ComponentContribution
from equilibrator_assets import (chemaxon, group_decompose, molecule,
                                 thermodynamics)
from equilibrator_cache import Compound, CompoundCache, CompoundMicrospecies
from openbabel.pybel import readstring
#KMS
from sqlalchemy import or_

logger = logging.getLogger(__name__)
group_decomposer = group_decompose.GroupDecomposer()
cc = ComponentContribution()
TRAINING_IDS = cc.predictor.params.train_G.index


def _populate_compound_information(row):
    """Attempt to populate a compound with key information.

    Accepts a pandas Series and attempts to generate a Compound.

    Returns a Compond or None if compound cannot be made.
    """
    try:

        if row.method == "chemaxon":
            # Use chemaxon to populate compound
            compound = Compound(**row.compound_dict)
            mol = molecule.Molecule.FromSmiles(compound.smiles)
        else:
            # Bypass chemaxon and use specified smiles
            compound_dictionary = {
                "atom_bag": chemaxon.get_atom_bag("smi", row.smiles),
                "dissociation_constants": [],
                "id": row.id,
                "smiles": row.smiles,
            }
            compound = Compound(**compound_dictionary)
            mol = molecule.Molecule.FromSmiles(row.smiles)

        if row.user_specified_pkas:
            # pkas in a string of form "A,B,C...""
            compound.dissociation_constants = [
                float(i) for i in row.user_specified_pkas.split(",")
            ]

        # Add extra information to compound
        compound.inchi_key = row.inchi_key
        try:
            decomposition = group_decomposer.Decompose(
                mol, ignore_protonations=False, raise_exception=True
            )
            compound.group_vector = decomposition.AsVector()
        except group_decompose.GroupDecompositionError:
            # Decomposition failed. If this is the first attempt
            # return None
            # If method is "empty" then store empty compound
            if row.method == "empty":
                compound.group_vector = None
            else:
                return None

        for ms_dict in thermodynamics.create_microspecies_mappings(compound):
            ms = CompoundMicrospecies(**ms_dict)
            compound.microspecies.append(ms)

        return compound
    
    except BaseException as e:
        print(e)
        return None


def create_compound(
    mol_strings: Union[str, List[str]],
    mol_format: str = "smiles",
    bypass_chemaxon: bool = False,
    save_empty_compounds: bool = False,
    specified_pkas: dict = None,
    log_df: pd.DataFrame = None,
    error_log: str = None,
) -> Union[Compound, List[Compound]]:
    """Generate a Compound object directly from SMILESs or InChIs.

    Parameters
    ----------
    mol_strings : Union[str, List[str]]
        Structure of compound(s) to add (InChI or smiles)
    mol_format : str, optional
        The format the molecules are given in (smiles, inchi),
        by default "smiles"
    bypass_chemaxon : bool, optional
        Allows compounds that fail to be decomposed with chemaxon to be
        created with the user-specified structure instead, by default False
    save_empty_compounds : bool, optional
            Whether or not to insert compounds into the database that cannot be
            decomposed user-specified structure, by default False
    specified_pkas : dict, optional
            A dictionary of user-specified pkas of form
            {mol_string: [pka1, pka2], ...}
            where mol_string is found in mol_strings, by default dict()
    log_df : pd.DataFrame, optional
            An empty pandas DataFrame to store the get_compounds results in to
            see the method compounds were obtained or if any fail,
            by default None
    error_log : str, optional
        File location to output any compounds that cannot be decomposed,
        by default None

    Returns
    -------
    Union[Compound, List[Compound]]
        The created compounds
    """
    if type(mol_strings) is str:
        _mol_strings = [mol_strings]
    else:
        _mol_strings = mol_strings

    molecules = pd.DataFrame(
        data=[[-1 - i, s] for i, s in enumerate(_mol_strings)],
        columns=[
            "id",
            "inchi",
        ],  # note that the "inchi" column can also contain SMILES strings
    )
    molecules["inchi_key"] = molecules.inchi.apply(
        lambda s: readstring(mol_format, s).write("inchikey").strip()
    )
    molecules["smiles"] = molecules.inchi.apply(
        lambda s: readstring(mol_format, s).write("smiles").strip()
    )
    molecules["compound_dict"] = list(
        thermodynamics.get_compound_mappings(
            molecules, "foo", num_acidic=20, num_basic=20
        )
    )
    # Specify dissociation constants if supplied
    molecules["user_specified_pkas"] = None
    if specified_pkas:
        for species, species_pkas in specified_pkas.items():
            if molecules["inchi"].str.contains(species, regex=False).any():
                molecules.loc[
                    molecules["inchi"] == species, ["user_specified_pkas"]
                ] = ",".join(map(str, species_pkas))

    # Find out how each compound can be succesfully inserted by sequentially
    # attempting to generate a compound with the following methods
    #   1. Using chemaxon to generate structure
    #   2. Using user-specified structure
    #   3. Inserting an empty compound
    molecules["compound"] = None
    for method in ["chemaxon", "bypass", "empty"]:
        molecules.loc[molecules["compound"].isnull(), "method"] = method
        molecules.loc[molecules["compound"].isnull(), "compound"] = molecules.loc[
            molecules["compound"].isnull()
        ].apply(_populate_compound_information, axis=1)

    # Rename inchi column, which isn't always inchi, to struct
    molecules = molecules.rename(columns={"inchi": "struct"})
    molecules["status"] = "valid"
    # Remove compounds from chemaxon and empty methods unless requested
    if not bypass_chemaxon:
        molecules.loc[molecules["method"] == "bypass", "compound"] = None
        molecules.loc[molecules["method"] == "bypass", "status"] = "failed"

    if not save_empty_compounds:
        molecules.loc[molecules["method"] == "empty", "compound"] = None
        molecules.loc[molecules["method"] == "empty", "status"] = "failed"

    # Log results
    molecules_string = molecules.to_string(
        columns=["struct", "inchi_key", "method", "status"]
    )
    if any(molecules["compound"].isnull()):
        logger.warning(
            "One or more compounds were unable to be decomposed."
            " Rerun specifying error_log or log_df to view details.\n"
        )
        logger.debug(f"Table of compound creation results\n{molecules_string}")
    else:
        logger.debug(
            "All compounds generated succesfully"
            f"Table of compound creation results\n{molecules_string}"
        )

    # Update log_df and output to error_log
    if (log_df is None) and error_log:
        log_df = pd.DataFrame()
        log_df[["struct", "inchi_key", "method", "status"]] = molecules[
            ["struct", "inchi_key", "method", "status"]
        ].values
    elif log_df is not None:
        common_ids = log_df.struct.isin(molecules.struct)
        log_df.loc[common_ids, ["inchi_key", "method", "status"]] = molecules[
            ["inchi_key", "method", "status"]
        ].values

    if error_log:
        log_df.to_csv(error_log, sep="\t", index=True)

    if type(mol_strings) is str:
        return molecules.compound.iat[0]
    else:
        return molecules.compound.tolist()


def get_or_create_compound(
    ccache: CompoundCache,
    mol_strings: Union[str, List[str]],
    mol_format: str = "smiles",
    connectivity_only: bool = False,
    bypass_chemaxon: bool = False,
    save_empty_compounds: bool = False,
    specified_pkas: dict = None,
    return_fails: bool = False,
    log_df: pd.DataFrame = None,
    error_log: str = None,
    read_only: bool = False,
) -> Union[Compound, List[Compound]]:
    """Get compounds from cache by descriptors, or creates them if missed.

    Parameters
    ----------
    ccache : CompoundCache
        [description]
    mol_strings : Union[str, List[str]]
        A string or list of strings containing text description of the
        molecule(s) (SMILES or InChI)
    mol_format : str, optional
        The format the molecules are given in (smiles, inchi),
        by default "smiles"
    connectivity_only : bool, optional
        Whether or not to use only connectivity portion of inchi_key when
        searching the ccache for an existing compound, by default False
    bypass_chemaxon : bool, optional
        Allows compounds that fail to be decomposed with chemaxon to be
        created with the user-specified structure instead, by default False
    save_empty_compounds : bool, optional
        Whether or not to insert compounds into the database that cannot be
        decomposed user-specified structure, by default False
    specified_pkas : dict, optional
            A dictionary of user-specified pkas of form
            {mol_string: [pka1, pka2], ...}
            where mol_string is found in mol_strings, by default dict()
    return_fails : bool, optional
        Whether or not to return failed compounds as None, by default False
    log_df : pd.DataFrame, optional
        An empty pandas DataFrame to store the get_compounds results in to
        see the method compounds were obtained or if any fail,
        by default None
    error_log : str, optional
        File location to output any compounds that cannot be decomposed,
        by default None
    read_only : bool
        Determines whether or not try attempt to create new compounds, or limit
        to existing compounds, by default False

    Returns
    -------
    Union[Compound, List[Compound]]
        Compound objects that were obtained from the database or created.
    """
    if type(mol_strings) is str:
        _mol_strings = [mol_strings]
    else:
        _mol_strings = mol_strings

    # create log_df if none was specified but error_log was
    if (log_df is None) and (error_log):
        log_df = pd.DataFrame()
    # allocate rows in dataframe with mol_strings
    if log_df is not None:
        log_df["struct"] = _mol_strings

    # InChI key is 3 parts separated by '-', X-Y-Z, where X is connectivity only
    # Y has stereochemical information, and Z describes deprotonation
    if connectivity_only:
        # Only take the connectivity block
        num_splits = 2
    else:
        # Take first two blocks
        num_splits = 1

    data = []
    for s in _mol_strings:
        inchi_key = readstring(mol_format, s).write("inchikey").strip()
        cc_search = ccache.search_compound_by_inchi_key(
            inchi_key.rsplit("-", num_splits)[0]
        )
        if cc_search:
            # Check if any compounds are in the training data
            training_compound = None
            for result in cc_search:
                if result.id in TRAINING_IDS:
                    training_compound = result
                    break

            if training_compound:
                # Found compound in training data
                data.append((s, training_compound))
            else:
                # No match, use lowest id number
                data.append((s, cc_search[0]))
        else:
            data.append((s, None))

    result_df = pd.DataFrame(data=data, columns=["mol_string", "compound"])
    misses = result_df.loc[result_df.compound.isnull(), :].index
    if not read_only:
        if len(misses) > 0:
            result_df.loc[misses, "compound"] = create_compound(
                result_df.loc[misses, "mol_string"].tolist(),
                mol_format,
                bypass_chemaxon,
                save_empty_compounds,
                specified_pkas,
                log_df,
            )
    else:
        result_df.loc[misses, "compound"] = None

    # Update the log
    if log_df is not None:
        success_ids = log_df.index.difference(misses)
        log_df.loc[success_ids, "method"] = "database"
        log_df.loc[success_ids, "status"] = "valid"
        log_df.loc[success_ids, "inchi_key"] = result_df.loc[
            success_ids, "compound"
        ].apply(lambda c: c.inchi_key)

        # Add in read-only status
        if read_only:
            log_df.loc[misses, "method"] = "read-only"
            log_df.loc[misses, "status"] = "failed"
            log_df.loc[misses, "inchi_key"] = None

    # Record log
    if error_log:
        log_df.to_csv(error_log, sep="\t", index=True)

    # Eliminate empty compounds from result_df if requested
    if not return_fails:
        result_df = result_df[~result_df["compound"].isnull()]

    if type(mol_strings) is str:
        # string as input -- return string or None
        return result_df.compound.iat[0] if not result_df.empty else None
    else:
        # list of strings as input -- return list
        return result_df.compound.tolist() if not result_df.empty else []

def get_or_create_compound_fast(
    ccache: CompoundCache,
    mol_strings: List[List[str]],
    mol_format: str = "smiles",
    connectivity_only: bool = False,
    bypass_chemaxon: bool = False,
    save_empty_compounds: bool = False,
    specified_pkas: dict = None,
    return_fails: bool = False,
    log_df: pd.DataFrame = None,
    error_log: str = None,
    read_only: bool = False,
) -> Union[Compound, List[Compound]]:
    """Get compounds from cache by descriptors, or creates them if missed.

    Parameters
    ----------
    ccache : CompoundCache
        [description]
    mol_strings : List[str, str]
        A list containing the SMILES and pickaxe_ids
        molecule(s) (SMILES or InChI)
    mol_format : str, optional
        The format the molecules are given in (smiles, inchi),
        by default "smiles"
    connectivity_only : bool, optional
        Whether or not to use only connectivity portion of inchi_key when
        searching the ccache for an existing compound, by default False
    bypass_chemaxon : bool, optional
        Allows compounds that fail to be decomposed with chemaxon to be
        created with the user-specified structure instead, by default False
    save_empty_compounds : bool, optional
        Whether or not to insert compounds into the database that cannot be
        decomposed user-specified structure, by default False
    specified_pkas : dict, optional
            A dictionary of user-specified pkas of form
            {mol_string: [pka1, pka2], ...}
            where mol_string is found in mol_strings, by default dict()
    return_fails : bool, optional
        Whether or not to return failed compounds as None, by default False
    log_df : pd.DataFrame, optional
        An empty pandas DataFrame to store the get_compounds results in to
        see the method compounds were obtained or if any fail,
        by default None
    error_log : str, optional
        File location to output any compounds that cannot be decomposed,
        by default None
    read_only : bool
        Determines whether or not try attempt to create new compounds, or limit
        to existing compounds, by default False

    Returns
    -------
    pd.DataFrame
        results_dataframe
    """
    # create log_df if none was specified but error_log was
    if (log_df is None) and (error_log):
        log_df = pd.DataFrame()
    # allocate rows in dataframe with mol_strings
    if log_df is not None:
        log_df["struct"] = mol_strings

    # InChI key is 3 parts separated by '-', X-Y-Z, where X is connectivity only
    # Y has stereochemical information, and Z describes deprotonation
    if connectivity_only:
        # Only take the connectivity block
        num_splits = 2
    else:
        # Take first two blocks
        num_splits = 1

    data = []
    
    mol_strings2pkid = {mol: pkid for mol, pkid in mol_strings}
    # Query things all at once instead of in for loop, the generate misses
    mol_strings2inchi = {mol: readstring(mol_format, mol).write("inchikey").strip().rsplit("-", num_splits)[0] for mol in mol_strings2pkid.keys()}
    mol_strings2inchi = {k: v for k, v in mol_strings2inchi.items() if v != ""}
    mols_inchi2strings = {v: k for k, v in mol_strings2inchi.items() if v != ""}
    inchi_queries = ccache.session.query(Compound).filter(
        or_(*[Compound.inchi_key.like(inchi_fragment + '%') for inchi_fragment in mol_strings2inchi.values()])
    ).all()

    inchi_dup_dict = defaultdict(list)
    for compound in inchi_queries:
        inchi_dup_dict[compound.inchi_key.rsplit("-", num_splits)[0]].append((compound.id, compound))

    mol_strings2compoundlist = {}
    # Sort the inchi_dup_list by id and then remove ids from list


    for inchi_dup_list in inchi_dup_dict.values():
        # Sort the list by id and select the lowest id
        sorted_compounds = [i[1] for i in sorted(inchi_dup_list, key=lambda x: x[0])]
        mol_strings2compoundlist[mols_inchi2strings[sorted_compounds[0].inchi_key.rsplit("-", num_splits)[0]]] = sorted_compounds
    for s, pkid in mol_strings:
        compound_list = mol_strings2compoundlist.get(s, None)
        if compound_list:
            training_compound = None
            for result in compound_list:
                if result.id in TRAINING_IDS:
                    training_compound = result
                    break

            if training_compound:
                # Found compound in training data
                data.append((s, training_compound, pkid))
            else:
                # No match, use lowest id number
                data.append((s, compound_list[0], pkid))      
        else:
            data.append((s, None, pkid))

    result_df = pd.DataFrame(data=data, columns=["mol_string", "compound", "pkid"])
    misses = result_df.loc[result_df.compound.isnull(), :].index
    if not read_only:
        if len(misses) > 0:
            result_df.loc[misses, "compound"] = create_compound(
                result_df.loc[misses, "mol_string"].tolist(),
                mol_format,
                bypass_chemaxon,
                save_empty_compounds,
                specified_pkas,
                log_df,
            )
    else:
        result_df.loc[misses, "compound"] = None

    # Update the log
    if log_df is not None:
        success_ids = log_df.index.difference(misses)
        log_df.loc[success_ids, "method"] = "database"
        log_df.loc[success_ids, "status"] = "valid"
        log_df.loc[success_ids, "inchi_key"] = result_df.loc[
            success_ids, "compound"
        ].apply(lambda c: c.inchi_key)

        # Add in read-only status
        if read_only:
            log_df.loc[misses, "method"] = "read-only"
            log_df.loc[misses, "status"] = "failed"
            log_df.loc[misses, "inchi_key"] = None

    # Record log
    if error_log:
        log_df.to_csv(error_log, sep="\t", index=True)

    # Eliminate empty compounds from result_df if requested
    if not return_fails:
        result_df = result_df[~result_df["compound"].isnull()]

    # list of strings as input -- return list
    return result_df if not result_df.empty else pd.DataFrame()
