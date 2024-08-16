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

'''
KMS CUSTOMIZATION OF EQUILIBRATOR
'''

from pathlib import Path
from typing import List, Union

import pandas as pd
from equilibrator_assets.chemaxon import get_chemaxon_status
from equilibrator_cache import Compound, CompoundIdentifier, Registry
from equilibrator_cache.api import create_compound_cache_from_sqlite_file
from equilibrator_cache.zenodo import (DEFAULT_COMPOUND_CACHE_SETTINGS,
                                       get_cached_filepath)

from .generate_compound import (get_or_create_compound,
                                get_or_create_compound_fast)

DEFAULT_CACHE_PATH = get_cached_filepath(DEFAULT_COMPOUND_CACHE_SETTINGS)

CHEMAXON_STATUS = get_chemaxon_status()

import time

from sqlalchemy import update


class LocalCompoundCache(object):
    """Read from and update a local compound cache."""

    def __init__(self, ccache_path: str = None):
        """Create a local cache object."""
        self.ccache = None
        self.ccache_path = None
        if ccache_path:
            self.load_cache(ccache_path)

        if CHEMAXON_STATUS == 0:  # cxcalc + license found
            self.read_only = False
            self._read_only_message = None
        elif CHEMAXON_STATUS == 1:  # cxcalc + no license
            self.read_only = True
            print(
                "No valid license for cxcalc installed, operating in read-only"
                " mode. A local cache may be loaded, but no compounds can be"
                " created. Please obtain a ChemAxon license to enable compound"
                " creation."
            )
        else:  # no cxcalc + no license
            self.read_only = True
            print(
                "cxcalc is not installed, operating in read-only"
                " mode. A local cache may be loaded, but no compounds can be"
                " created. Install cxcalc and obtain a ChemAxon license to"
                " enable compound creation."
            )

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
        elif ccache_path == DEFAULT_COMPOUND_CACHE_SETTINGS.filename:
            print(
                "Default eQuilibrator cache cannot be used with"
                "LocalCompoundCache. Make a cache copy using"
                "load_compound_cache_from_zenodo."
            )
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

        elif new_ccache_path == DEFAULT_CACHE_PATH:
            print(
                "Default eQuilibrator compound cache specified." " Specify a new file."
            )
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

        print(f"Copying default Zenodo compound cache to {new_ccache_path}")
        new_ccache_path.write_bytes(DEFAULT_CACHE_PATH.read_bytes())

    def get_compounds(
        self,
        mol_strings: Union[str, List[str]],
        mol_format: str = "smiles",
        connectivity_only: bool = False,
        bypass_chemaxon: bool = False,
        save_empty_compounds: bool = False,
        specified_pkas: dict = None,
        return_fails: bool = False,
        log_df: pd.DataFrame = None,
        error_log: str = None,
    ) -> Union[Compound, List[Compound]]:
        """Get the Compound object of a list of molecules.

        Get compounds from the CompoundCache. If any compounds are not found,
        attempt to create them and insert them into the CompoundCache.

        Parameters
        ----------
        mol_strings : Union[str, List[str]]
            Structure of compound(s) to add (InChI or smiles)
        mol_format : str, optional
            The format the molecules are given in ("inchi" or "smiles"),
            by default "smiles"
        connectivity_only : bool, optional
            Whether or not to use the connectivity only portion of the
            inchi-key to search, by default False
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
            Whether or not to omit or return failed compounds as None,
            by default False
        log_df : pd.DataFrame, optional
            An empty pandas DataFrame to store the get_compounds results in to
            see the method compounds were obtained or if any fail,
            by default None
        error_log : str, optional
            File location (.tsv) to write any compounds that cannot be
            decomposed, by default None

        Returns
        -------
        Union[Compound, List[Compound]]
            Compound objects that were obtained from the database or created.
        """
        if not self.ccache:
            print("No cache found: load a cache with load_cache() first.")
            return None

        if self.read_only:
            print(
                "Read-Only mode: Only existing compounds"
                " in the database can be accessed."
            )

        synonym_registry = (
            self.ccache.session.query(Registry).filter_by(namespace="synonyms").one()
        )

        # create log_df if none was specified but error_log was
        if (log_df is None) and (error_log):
            log_df = pd.DataFrame()

        # print("get_or_create_compound_start", flush=True)
        start1 = time.time()
        compounds = get_or_create_compound(
            self.ccache,
            mol_strings,
            mol_format,
            connectivity_only,
            bypass_chemaxon,
            save_empty_compounds,
            specified_pkas,
            return_fails,
            log_df,
            error_log,
            self.read_only,
        )
        # print("get_or_create_compound_end", time.time() - start1, flush=True)

        # Record log
        if error_log:
            log_df.to_csv(error_log, sep="\t", index=True)

        # Process compound results.
        # If compounds have a negative id do the following
        # for insertion to .sqlite:
        #  1. Delete ID so insertion to db assigns automatic ID
        #  2. change group_vector to a list to avoid pickle issues upon recall
        #  3. Add a default synonym that matches compound.id
        if type(compounds) is not list:
            compounds = [compounds]

        # print("start_adding_compounds", flush=True)
        start1 = time.time()
        # Add in new compounds and keep track of indices
        new_compounds = []
        new_objects = []
        for i, compound in enumerate(compounds):
            if compound:
                if compound.id <= -1:
                    new_compounds.append(i)
                    del compound.id
                    # Make group_vec list for de-pickeling
                    if compound.group_vector:
                        compound.group_vector = list(compound.group_vector)
                    # insert compound to local 
                    new_objects.append(compound)
        
        # insert compounds to local session
        self.ccache.session.bulk_save_objects(new_objects)

        # Commit to get automatically generated ID
        self.ccache.session.commit()
        # print("add_compounds_end", time.time() - start1, flush=True)

        # Add a default synonym for new compounds
        # print("start_adding_synonyms", flush=True)
        start1 = time.time()
        new_objects = []
        for idx in new_compounds:
            compound = compounds[idx]
            new_objects.append(
                CompoundIdentifier(
                    registry=synonym_registry,
                    accession=compound.id,
                    compound_id=compound.id,
                )
            )

        # insert identifiers to local session
        self.ccache.session.bulk_save_objects(new_objects)
        self.ccache.session.commit()
        # print("add_synonyms_end", time.time() - start1, flush=True)

        if type(mol_strings) is str:
            compounds = compounds[0]

        return compounds
    

    def get_compounds_fast(
            self,
            mol_strings: List[List[str]],
            mol_format: str = "smiles",
            connectivity_only: bool = False,
            bypass_chemaxon: bool = False,
            save_empty_compounds: bool = False,
            specified_pkas: dict = None,
            return_fails: bool = False,
            log_df: pd.DataFrame = None,
            error_log: str = None,
        ) -> Union[Compound, List[Compound]]:
            """Get the Compound object of a list of molecules.

            Get compounds from the CompoundCache. If any compounds are not found,
            attempt to create them and insert them into the CompoundCache.

            Parameters
            ----------
            mol_strings : List[str]
                A list of SMILES and pickaxe IDs
            mol_format : str, optional
                The format the molecules are given in ("inchi" or "smiles"),
                by default "smiles"
            connectivity_only : bool, optional
                Whether or not to use the connectivity only portion of the
                inchi-key to search, by default False
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
                Whether or not to omit or return failed compounds as None,
                by default False
            log_df : pd.DataFrame, optional
                An empty pandas DataFrame to store the get_compounds results in to
                see the method compounds were obtained or if any fail,
                by default None
            error_log : str, optional
                File location (.tsv) to write any compounds that cannot be
                decomposed, by default None

            Returns
            -------
            Union[Compound, List[Compound]]
                Compound objects that were obtained from the database or created.
            """
            if not self.ccache:
                print("No cache found: load a cache with load_cache() first.")
                return None

            if self.read_only:
                print(
                    "Read-Only mode: Only existing compounds"
                    " in the database can be accessed."
                )

            coco_registry = (
                self.ccache.session.query(Registry).filter_by(namespace="coco").one()
            )

            # create log_df if none was specified but error_log was
            if (log_df is None) and (error_log):
                log_df = pd.DataFrame()

            import time
            start1 = time.time()
            # print("get_or_create_compound_fast_start", flush=True)
            compounds_df = get_or_create_compound_fast(
                self.ccache,
                mol_strings,
                mol_format,
                connectivity_only,
                bypass_chemaxon,
                save_empty_compounds,
                specified_pkas,
                return_fails,
                log_df,
                error_log,
                self.read_only,
            )
            # print("get_or_create_compound_fast", time.time() - start1, flush=True)

            # Record log
            if error_log:
                log_df.to_csv(error_log, sep="\t", index=True)

            # Process compound results.
            # If compounds have a negative id do the following
            # for insertion to .sqlite:
            #  1. Delete ID so insertion to db assigns automatic ID
            #  2. change group_vector to a list to avoid pickle issues upon recall
            #  3. Add a default synonym that matches compound.id


            # print("start_adding_compounds", flush=True)
            start2 = time.time()
            # Add in new compounds and keep track of indices
            new_compounds = []
            new_objects = []
            for i, row in compounds_df.iterrows():
                if row.compound:
                    if row.compound.id <= -1:
                        new_compounds.append(i)
                        del(row.compound.id)
                        # Make group_vec list for de-pickeling
                        if row.compound.group_vector:
                            row.compound.group_vector = list(row.compound.group_vector)
                
                # print("added identifier", time.time() - start2, flush=True)

                # print("adding individual, adding new objects")
                new_objects.append(row.compound)
                # print("added_new_objects", time.time() - start2, flush=True)

            # insert compounds to local session
            self.ccache.session.add_all(new_objects)
            self.ccache.session.commit()
            # print("add_compounds", time.time() - start2, flush=True)
            
            # print("start_adding_coco_registries", flush=True)
            new_objects = []
            for i, row in compounds_df.iterrows():
                if row.compound:
                    # print("adding individual, registering identifier", time.time() - start2, flush=True)     
                    compound_identifier = CompoundIdentifier(
                            accession=row.pkid,
                            compound_id=row.compound.id,
                        )     
                    new_objects.append(compound_identifier)

            # insert identifiers to local session
            self.ccache.session.bulk_save_objects(new_objects)
            self.ccache.session.commit()
            # print("added_coco_registries", time.time() - start2, flush=True)
            
            # print("start_updating_coco_registries", flush=True)
            stmt = (
                update(CompoundIdentifier).
                where(CompoundIdentifier.accession.in_([row.pkid for _, row in compounds_df.iterrows()])).
                values(registry_id=coco_registry.id)
            )
            self.ccache.session.execute(stmt)
            self.ccache.session.commit()

            # print("added_coco_registries", time.time() - start2, flush=True)      
            # Commit to get automatically generated ID
            # print("commit", flush=True)
            startc = time.time()
            self.ccache.session.commit()
            # print("commitend", time.time() - startc, flush=True)
            
            
            # No need for pickaxe, regular does need this
            # print("start_adding_synonyms", flush=True)
            # start3 = time.time()
            # for idx in new_compounds:
            #     compound = compounds_df.iloc[idx].compound
            #     compound.identifiers.append(
            #         CompoundIdentifier(
            #             registry=synonym_registry,
            #             accession=compound.id,
            #             compound_id=compound.id,
            #         )
            #     )

            # print("add_synonyms", time.time() - start3, flush=True)

            startc = time.time()
            # print("commit", flush=True)
            self.ccache.session.commit()
            # print("commitend", time.time() - startc, flush=True)

            if type(mol_strings) is str:
                compounds = compounds[0]

            return compounds_df.compound.tolist()
    
    def add_pickaxe_compounds(
        self,
        compound_df: pd.DataFrame,
        mol_format: str = "smiles",
        connectivity_only: bool = False,
        bypass_chemaxon: bool = False,
        save_empty_compounds: bool = False,
        specified_pkas: dict = None,
        log_df: pd.DataFrame = None,
        error_log: str = None,
    ) -> None:
        """Add a dataframe of compounds and ids to the cache for pickaxe specifically.

        This code will use a faster query to check if a compound exists in the database

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
            by default "compound_creation_log.tsv"
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

        # Attempt to add every compound to the local cache
        # Afterwards insert coco_id if not already existing
        #TODO remove
        import time
        start1 = time.time()
        # print("get_compounds_fast_start", flush=True)
        compound_df["compound"] = self.get_compounds_fast(
            compound_df[["struct", "name"]].values.tolist(),
            mol_format,
            connectivity_only,
            bypass_chemaxon,
            save_empty_compounds,
            specified_pkas,
            True,  # Return fails
            log_df,
            error_log,
        )
        # print("get_compounds_fast_end", time.time() - start1, flush=True)

        # print("start_adding_compounds", flush=True)
        start2 = time.time()

        # Create Compound Identifiers
        # for row in compound_df.itertuples(index=False):
        #     print("starting_new_cpd", time.time()-start2, flush=True)
        #     if row.compound:
        #         coco_identifier = None
        #         synonym = None

        #         if row.coco_id:
        #             coco_identifier = row.coco_id
        #         elif row.name:
        #             coco_identifier = row.name
        #         # print("starting_single_coco", time.time()-start2, flush=True)
        #         # if coco_identifier:
        #         #     coco_identifier = CompoundIdentifier(
        #         #         registry=coco_registry,
        #         #         accession=coco_identifier,
        #         #         compound_id=row.compound.id,
        #         #     )
        #         #     # In regular we check to see if it already exists, for pickaxe don't
        #         #     row.compound.identifiers.append(coco_identifier)
        #         #     print("added_single_coco", time.time()-start2, flush=True)

        #         # Assign synonym identifier, use coco_id if name isn't
        #         # available and finally use compound ID if neither name
        #         # nor coco_id is available.
        #         if row.name:
        #             synonym = row.name
        #         elif row.coco_id:
        #             synonym = row.coco_id

        #         if synonym:
        #             remove_default_synonym(row.compound)
        #             synonym = CompoundIdentifier(
        #                 registry=synonym_registry,
        #                 accession=synonym,
        #                 compound_id=row.compound.id,
        #             )
        #             # In regular we check to see if it already exists, for pickaxe don't
        #             row.compound.identifiers.append(synonym)
        #             print("added_single_synonym", time.time()-start2, flush=True)
            
            

        # self.ccache.session.commit()

        # print("add_compounds", time.time() - start2, flush=True)

    def add_compounds(
        self,
        compound_df: pd.DataFrame,
        mol_format: str = "smiles",
        connectivity_only: bool = False,
        bypass_chemaxon: bool = False,
        save_empty_compounds: bool = False,
        specified_pkas: dict = None,
        log_df: pd.DataFrame = None,
        error_log: str = None,
    ) -> None:
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
            by default "compound_creation_log.tsv"
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

        # Attempt to add every compound to the local cache
        # Afterwards insert coco_id if not already existing
        # print("get_compounds_start", flush=True)
        start1 = time.time()
        compound_df["compound"] = self.get_compounds(
            list(compound_df["struct"]),
            mol_format,
            connectivity_only,
            bypass_chemaxon,
            save_empty_compounds,
            specified_pkas,
            True,  # Return fails
            log_df,
            error_log,
        )
        # print("get_compounds_end", time.time() - start1, flush=True)

        # print("start_adding_compounds", flush=True)
        start2 = time.time()
        # Create Compound Identifiers
        new_objects = []
        for row in compound_df.itertuples(index=False):
            if row.compound:
                # print("starting_single_compound", time.time()-start2, flush=True)
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
                        compound_id=row.compound.id,
                    )
                    if not in_identifiers(row.compound, coco_identifier):
                        new_objects.append(coco_identifier)
                        # row.compound.identifiers.append(coco_identifier)

                    # print("added_single_coco", time.time()-start2, flush=True)
                # Assign synonym identifier, use coco_id if name isn't
                # available and finally use compound ID if neither name
                # nor coco_id is available.
                if row.name:
                    synonym = row.name
                elif row.coco_id:
                    synonym = row.coco_id

                if synonym:
                    remove_default_synonym(row.compound)
                    synonym = CompoundIdentifier(
                        registry=synonym_registry,
                        accession=synonym,
                        compound_id=row.compound.id,
                    )
                    if not in_identifiers(row.compound, synonym):
                        new_objects.append(synonym)
                    # print("added_single_synonym", time.time()-start2, flush=True)
                # print("added_single_compound", time.time()-start2, flush=True)

            # print("starting_bulk", time.time()-start2, flush=True)
            start3 = time.time()
            self.ccache.session.bulk_save_objects(new_objects)
            self.ccache.session.commit()
            # print("bulk end", time.time()-start3, flush=True)

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
