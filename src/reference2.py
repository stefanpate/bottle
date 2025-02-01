from typing import Callable, List
from minedatabaseV2.filters.base_filter import Filter
from minedatabaseV2.pickaxe import Pickaxe
from minedatabaseV2.filters.base_filter import Filter
import requests
import copy
import multiprocessing
import time
from functools import partial
from typing import Callable, List, Set, Tuple
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from rdkit import Chem
from minedatabaseV2.utils import Chunks
logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

class KnownFilter(Filter):
    """Filter that checks if intermediate metabolites are present within the biological databases BRENDA, KEGG, METACYC
    ----------

    Attributes
    ----------
    """

    def __init__(self, filepath:str) -> None:
        self._filter_name = "Known Biological Compounds Filter"
        self.filepath = filepath
        self.known_cpds_db = set(line.strip() for line in open(self.filepath))

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _pre_print(self) -> None:
        """Print before filtering."""
        print(
            (f"Searching for known compounds within biological databasese"))

    def _post_print(self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float) -> None:
        """Print after filtering."""
        print((f"{n_filtered} of {n_total} "
                "compounds selected after "
                f"Similarity Sampling of generation {pickaxe.generation}"
                f"--took {time.time() - time_sample}s.\n"))

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """
        Samples N compounds to expand based on the weighted Similarity distribution.

        Parameters
        ----------
        pickaxe : Pickaxe
            Pickaxe object to filter
        processes : int
            Number of processes to use.
        """

        def canonicalize_smiles(smi: str):
            """
            Canonicalize SMILES string with RDKit
            :param: Input SMILES string, assume not canonical
            :return: Canonicalized SMILES string
            """

            uncanon_smi = smi  # assume input SMILES string is not canonical

            try:  # try canonicalizing first
                canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(uncanon_smi))

            except:  # leave SMILE as is if it cannot be canonicalized
                canon_smi = uncanon_smi

            return canon_smi

        def search_SMILES_in_biological_databases(SMILES: str, biological_SMILES: set):
            """
            Search SMILES of a compound against combined BRENDA, KEGG, and Metacyc databases
            @param SMILES: SMILES of queried compound
            @param biological_SMILES: list of SMILES present in BRENDA, KEGG, and Metacyc databases
            @return: True if compound is present in biological databases
            """
            canonicalized_SMILES = canonicalize_smiles(SMILES)  # canonicalize queried SMILES
            return (
                    canonicalized_SMILES in biological_SMILES
            )  # True if compound is present in biological databases

        def smiles_to_pubchem_cid(smiles):
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                return data['IdentifierList']['CID'][0]
            else:
                return None

        def is_target(cpd, pickaxe):
            for t_id in pickaxe.targets:
                if "C" + t_id[1:] == cpd["_id"]:
                    return True
            return False

        print(f"Filtering Generation {pickaxe.generation} by searching for known compounds within Pubchem")

        cpds_remove_set = set()
        rxn_remove_set = set()

        # Get compounds eligible for expansion in the current generation
        compounds_to_check = []
        set_unreactive = False

        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation \
                    and pickaxe.generation != 0 \
                    and cpd["Type"] not in ["Coreactant","Target Compound"]:

                if not smiles_to_pubchem_cid(cpd["SMILES"]):
                    cpds_remove_set.add(cpd["_id"])
                    pickaxe.compounds[cpd["_id"]]["Expand"] = False

                # if not search_SMILES_in_biological_databases(cpd["SMILES"], self.known_cpds_db):
                #     cpds_remove_set.add(cpd["_id"])
                #     pickaxe.compounds[cpd["_id"]]["Expand"] = False

                # #print(search_SMILES_in_biological_databases(cpd["SMILES"],self.known_cpds_db))
                # # Check for targets and only react if terminal
                # print(cpd,pickaxe.react_targets)
                # if pickaxe.react_targets:
                #     compounds_to_check.append(cpd)
                # else:
                #     print(cpd)
                #     if is_target(cpd, pickaxe):
                #         pickaxe.compounds[cpd["_id"]]["Expand"] = False
                #     else:
                #         if not search_SMILES_in_biological_databases(cpd["SMILES"], self.known_cpds_db):
                #             cpds_remove_set.add(cpd["_id"])

        for c_id in cpds_remove_set:
            pickaxe.compounds[c_id]["Expand"] = False

        return cpds_remove_set, rxn_remove_set