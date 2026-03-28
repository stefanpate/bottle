"""Store and update a local compound cache."""

# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
# Copyright (c) 2018 Institute for Molecular Systems Biology,
# ETH Zurich, Switzerland.
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

import logging
from pathlib import Path
from typing import List

import pandas as pd
from equilibrator_cache import Compound, CompoundIdentifier, Registry
from equilibrator_cache.api import create_compound_cache_from_sqlite_file
from equilibrator_assets.generate_compound import GenerateCompoundResult, get_or_create_compounds


logger = logging.getLogger(__name__)


class LocalCompoundCache(object):
    """Read from and update a local compound cache."""

    def __init__(self, ccache_path: str = None):
        """Create a local cache object."""
        self.ccache = None
        self.ccache_path = None
        if ccache_path:
            self.load_cache(ccache_path)

    def load_cache(self, ccache_path: str) -> None:
        """Load a cache from a .sqlite file locally.

        Load a local cache from a compound.sqlite file derived from
        the eQuilibrator Zenodo data.

        Parameters
        ----------
        ccache_path : str
            The location from which to load the cache.
        """
        ccache_path = Path(ccache_path)

        if ccache_path.suffix != ".sqlite":
            print("Provided file is not a .sqlite file.")
        elif ccache_path.is_file():
            self.ccache_path = ccache_path
            print(f"Loading compounds from {ccache_path}")
            if self.ccache:
                self.ccache.session.close()
            self.ccache = create_compound_cache_from_sqlite_file(ccache_path)
        else:
            print(f"{ccache_path} does not exist.")

        return None

    @staticmethod
    def generate_local_cache_from_default_zenodo(
        new_ccache_path: str, force_write: bool = False
    ) -> None:
        """Create a local cache from the default zenodo .sqlite.

        Parameters
        ----------
        new_ccache_path : str
            The folder to export the cache into.

        force_write : bool
            Write new cache by overwriting and generating a new directory
            if necessary, default False
        """
        new_ccache_path = Path(new_ccache_path)

        # Ensure file type and don't allow overwriting of default cache
        if new_ccache_path.suffix != ".sqlite":
            print("New compound cache must be a .sqlite extension.")
            return None

        if force_write:
            # Make sure folder exists and file is deleted for writing
            if not new_ccache_path.parent.is_dir():
                new_ccache_path.parent.mkdir()
            elif new_ccache_path.is_file():
                new_ccache_path.unlink()

        elif new_ccache_path.is_file():
            print(f"{new_ccache_path} already exists.")
            print("Delete existing file and replace?")

            choice = ""
            while choice.lower() not in ["yes", "no"]:
                choice = input("Proceed? (yes/no):")

            if choice.lower() == "no":
                print("Local cache generation cancelled.")
                return None
            else:
                print(f"Deleting {new_ccache_path}")
                new_ccache_path.unlink()

        elif not new_ccache_path.parent.is_dir():
            print(f"{new_ccache_path.parent} does not exist. Create?")
            choice = ""
            while choice.lower() not in ["yes", "no"]:
                choice = input("Proceed? (yes/no):")

            if choice.lower() == "yes":
                new_ccache_path.parent.mkdir()
            else:
                return None

    def get_compounds(
        self,
        mol_strings: List[str],
        mol_format: str = "smiles",
        connectivity_only: bool = False,
        bypass_chemaxon: bool = False,
        save_empty_compounds: bool = False,
        specified_pkas: dict = None,
    ) -> List[GenerateCompoundResult]:
        """Get the Compound object of a list of molecules.

        Get compounds from the CompoundCache. If any compounds are not found,
        attempt to create them and insert them into the CompoundCache.

        Parameters
        ----------
        mol_strings : List[str]
            Structure of compounds to add (InChI or smiles)
        mol_format : str, optional
            The format the molecules are given in ("inchi" or "smiles"),
            by default "smiles"
        connectivity_only : bool, optional
            Whether to use the connectivity only portion of the
            inchi-key to search, by default False
        bypass_chemaxon : bool, optional
            Allows compounds that fail to be decomposed with chemaxon to be
            created with the user-specified structure instead, by default False
        save_empty_compounds : bool, optional
            Whether to insert compounds into the database that cannot be
            decomposed user-specified structure, by default False
        specified_pkas : dict, optional
            A dictionary of user-specified pkas of form
            {mol_string: [pka1, pka2], ...}
            where mol_string is found in mol_strings, by default dict()

        Returns
        -------
        List[GenerateCompoundResult]
            Compound objects that were obtained from the database or created.
        """
        if not self.ccache:
            logger.debug("No cache found: load a cache with load_cache() first.")
            return None

        if self.read_only:
            logger.debug(
                "Read-Only mode: Only existing compounds"
                " in the database can be accessed."
            )

        synonym_registry = (
            self.ccache.session.query(Registry).filter_by(namespace="synonyms").one()
        )

        cpd_results = get_or_create_compounds(
            ccache=self.ccache,
            mol_strings=mol_strings,
            mol_format=mol_format,
            connectivity_only=connectivity_only,
            bypass_chemaxon=bypass_chemaxon,
            save_empty_compounds=save_empty_compounds,
            specified_pkas=specified_pkas,
            read_only=self.read_only,
        )

        # Add in new compounds and keep track of indices
        new_compounds = []
        for i, cpd_res in enumerate(cpd_results):
            compound = cpd_res.compound
            if compound:
                if compound.id <= -1:
                    new_compounds.append(i)
                    del compound.id
                    # Make group_vec list for de-pickeling
                    if compound.group_vector:
                        compound.group_vector = list(compound.group_vector)
                    # insert compound to local session
                    self.ccache.session.add(compound)

        # Commit to get automatically generated ID
        self.ccache.session.commit()

        # Add a default synonym for new compounds
        for idx in new_compounds:
            compound = cpd_results[idx].compound
            compound.identifiers.append(
                CompoundIdentifier(
                    registry=synonym_registry,
                    accession=compound.id,
                    compound_id=compound.id,
                )
            )

        self.ccache.session.commit()

        return cpd_results

    def add_compounds(
        self,
        compound_df: pd.DataFrame,
        mol_format: str = "smiles",
        connectivity_only: bool = False,
        bypass_chemaxon: bool = False,
        save_empty_compounds: bool = False,
        specified_pkas: dict = None,
    ) -> List[GenerateCompoundResult]:
        """Add a dataframe of compounds and ids to the cache.

        Takes a dataframe with three columns:
            struct -- the compound structures in either smiles or inchi
            coco_id -- the id of the compounds to insert into the coco namespace
            name -- a common name for the compound

        An attempt is then made to generate the compound and insert it into the
        database. Depending on the input values the function will add specific
        values into the database.

        1. struct only: The compound is assigned a default ID that is used as
        its ID in the coco namespace and as its common name.
        2. struct and coco_id: The compound is assigned the coco_id as its
        coco identifier and common name.
        3. struct and name: The compound is assigned the name as its
        coco identifier and common name.
        4. struct, id, and name: The compound is assigned the coco_id as its
        coco identifier and assigned thename as the common name.

        Parameters
        ----------
        compound_df : pd.DataFrame
            A pandas dataframe with two columns, "struct", "id", and "name"
        mol_format : str, optional
            The format the molecules are given in ("smiles" or "inchi"),
            by default "smiles"
        connectivity_only : bool, optional
            Whether or not to use the connectivity only portion of the
            inchi-key to search, by default False
        bypass_chemaxon : bool, optional
            Allows compounds that fail to be decomposed with chemaxon to be
            created with the user-specified structure instead, by default False
        save_empty_compounds : bool, optional
            Whether to insert compounds into the database that cannot be
            decomposed user-specified structure, by default False
        specified_pkas : dict, optional
            A dictionary of user-specified pkas of form
            {mol_string: [pka1, pka2], ...}
            where mol_string is found in mol_strings, by default dict()
        """
        # Check if specified coco_id is already in compounds registry
        if not self.ccache:
            print("No cache found: load a cache with load_cache() first.")
            return None

        def in_identifiers(compound, new_identifier):
            for identifier in compound.identifiers:
                if (
                    identifier.registry == new_identifier.registry
                    and identifier.accession == new_identifier.accession
                ):
                    return True
                else:
                    continue
            return False

        def remove_default_synonym(compound):
            # The default synonym matches the compound id
            # Check for matches and deleted
            for identifier in compound.identifiers:
                if (
                    identifier.accession.isdigit()
                    and identifier.registry == synonym_registry
                ):
                    if int(identifier.accession) == compound.id:
                        self.ccache.session.delete(identifier)
                        return None

        coco_registry = (
            self.ccache.session.query(Registry).filter_by(namespace="coco").one()
        )

        synonym_registry = (
            self.ccache.session.query(Registry).filter_by(namespace="synonyms").one()
        )

        # Attempt to add every compound to the local cache.
        # Afterwards insert coco_id if not already existing
        cpd_results = self.get_compounds(
            mol_strings=compound_df["struct"].to_list(),
            mol_format=mol_format,
            connectivity_only=connectivity_only,
            bypass_chemaxon=bypass_chemaxon,
            save_empty_compounds=save_empty_compounds,
            specified_pkas=specified_pkas,
        )

        # Create Compound Identifiers
        for i, row in enumerate(compound_df.itertuples()):
            compound = cpd_results[i].compound
            if compound is not None:
                coco_identifier = None
                synonym = None

                # Create coco_id identifier if possible
                # Prioritize specified coco_id, but try name if coco_id
                # is unavailable
                if row.coco_id:
                    coco_identifier = row.coco_id
                elif row.name:
                    coco_identifier = row.name

                if coco_identifier:
                    coco_identifier = CompoundIdentifier(
                        registry=coco_registry,
                        accession=coco_identifier,
                        compound_id=compound.id,
                    )
                    if not in_identifiers(compound, coco_identifier):
                        compound.identifiers.append(coco_identifier)

                # Assign synonym identifier, use coco_id if name isn't
                # available and finally use compound ID if neither name
                # nor coco_id is available.
                if row.name:
                    synonym = row.name
                elif row.coco_id:
                    synonym = row.coco_id

                if synonym:
                    remove_default_synonym(compound)
                    synonym = CompoundIdentifier(
                        registry=synonym_registry,
                        accession=synonym,
                        compound_id=compound.id,
                    )
                    if not in_identifiers(compound, synonym):
                        compound.identifiers.append(synonym)

            self.ccache.session.commit()

        return cpd_results

    def get_coco_accessions(self) -> List[str]:
        """Return all accessions in the coco namespace.

        Returns
        -------
        List[str]
            A list of accessions for compounds in the coco namespace.
        """
        query = (
            self.ccache.session.query(Compound)
            .outerjoin(CompoundIdentifier)
            .outerjoin(Registry)
        )

        compounds = query.filter(Registry.namespace == "coco").all()

        coco_ids = []
        for compound in compounds:
            for compound_id in compound.identifiers:
                if compound_id.registry.namespace == "coco":
                    coco_ids.append(compound_id.accession)

        return coco_ids
