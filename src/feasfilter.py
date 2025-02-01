import abc
import time
from copy import copy
from typing import List, Set
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from typing import Callable, List
import os 
import pickle
import numpy as np 
import pandas as pd
from DORA_XGB import DORA_XGB
import json

logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

from DORA_XGB import DORA_XGB
from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe

######################################
 # DORAXGB

class DORAXGBFilter(Filter):

    def __init__(self, method:str, generation_list =[],last_generation_only=False,threshold:str,) -> None:
        self._filter_name = "DORAXGB Reaction Filter"
        self.method = method 
        self.threshold = threshold
        self.generation_list = generation_list
        self.last_generation_only = False


    def filter_name(self)-> str:
        return self._filter_name
    
    def _pre_print(self) -> None:
        """Print filter being applied."""
        print(f"Applying filter: {self.filter_name}")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float) -> None:
        """Print results of filtering.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
            Unused here, but may be useful in your implementation.
        n_total : int
            Total number of compounds.
        n_filtered : int
            Number of compounds remaining after filtering.
        times_sample : float
            Time in seconds from time.time().
        """
        print((f"{n_filtered} of {n_total}"
            "reactions remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n")
        )
    
    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Set[str]:
        """Return list of reactions to remove from pickaxe object.

        Parameters
        ----------
        pickaxe : Pickaxe
            Instance of Pickaxe being used to expand and filter the network.
        processes : int
            The number of processes to use, by default 1.
        generation : int
            Which generation the expansion is in.
        """
        rxn_remove_set = set()
    
        def feasibility_score(rxn_id: str) -> float:
            """Returns the feasibility score for a reaction. """
            by_desc_MW_model = DORA_XGB.feasibility_classifier(cofactor_positioning = 'by_descending_MW')
            products = pickaxe.reactions[rxn_id]["Products"]
            score = by_desc_MW_model.predict_proba(products)  # Get feasibility score
            return score  # Simply return the score
        
        # def should_delete_reaction(rxn_ids: str, threshold: float) -> bool:
        #     """
        #     Returns True if the reaction should be removed (if any product's score is below threshold).
        #     Returns False if the reaction is valid (all product scores are above threshold).
        #     """
        #     scores = feasibility_score(rxn_id)  # Get product-wise scores for the reaction
        #     # Check if any product's score is below the threshold
        #     if any (score < threshold for score in scores.values()):
        #         return True 
        #     else:
        #         return False

        # def remove_reaction(pickaxe, rxn_ids: list, threshold: float) -> set:
        #     """
        #     Removes reactions where any product's score is below the threshold.
    
        #     Args:
        #     pickaxe: The Pickaxe object containing reactions.
        #     rxn_ids: List of reaction IDs to evaluate.
        #     threshold: The threshold score below which reactions should be removed.
    
        #     Returns:
        #     A set of reaction IDs that passed the feasibility check (valid reactions).
        #     """
        #     reactions_to_return = set()

        #     for rxn_id in list(rxn_ids):  # Use list() to avoid modifying the dictionary while iterating
        #         products = pickaxe.reactions[rxn_id]["Products"]
        #         scores = feasibility_score(rxn_id)  # Get product-wise scores for the reaction
        #         # Check if all product scores meet or exceed the threshold
        #         if all(score >= threshold for score in scores.values()):
        #             reactions_to_return.add(rxn_id)
        #         else:
        #             print(f"Reaction {rxn_id} removed due to low feasibility score.")
        #             pickaxe.reactions.pop(rxn_id, None)  # Safe deletion
        #             #Store removed reaction details
        #             reaction_ids_to_delete()

        #     return reactions_to_return
        
if __name__ == "__main__":
    pass