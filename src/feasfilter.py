
# --------------------------------- WORKS keeping the compounds "Expand" = True------------------------------------
import time
import logging
import numpy as np
from typing import Set, Tuple
import rdkit.rdBase as rkrb
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
import rdkit.RDLogger as rkl
from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe
from DORA_XGB import DORA_XGB

# Suppress RDKit warnings
logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

class DORAXGBFilter(Filter):
    """
    A filter that removes reactions and associated compounds based on a feasibility score.
    
    This filter uses a machine learning model to assess the feasibility of reactions 
    and removes those that do not meet the specified feasibility threshold. 
    It also removes compounds associated with the removed reactions unless they are 
    of a special type (Coreactant, Target Compounds, Starting Compound).
    
    """
    
    def __init__(self, threshold: float, generation_list=None, last_generation_only=False) -> None:
        """
        Initializes the DORAXGBFilter with the specified parameters.

        Args:
            threshold (float): The threshold below which reactions are considered unfeasible.
            generation_list (list, optional): A list of generations to filter (default is all generations).
            last_generation_only (bool, optional): If True, only the last generation will be filtered (default is False).
        """
        self._filter_name = "DORAXGB Reaction Filter"
        self.threshold = threshold
        self.generation_list = generation_list or []
        self.last_generation_only = last_generation_only
        self.model = DORA_XGB.feasibility_classifier(cofactor_positioning="add_concat", model_type="spare")

    def filter_name(self) -> str:
        """
        Returns the name of the filter.

        Returns:
            str: The name of the filter.
        """
        return self._filter_name

    def _pre_print_header(self, pickaxe: Pickaxe) -> None:
        """
        Prints the header before starting the filtering process.

        Args:
            pickaxe (Pickaxe): The object representing the current pickaxe with reactions and compounds.
        """
        print("----------------------------------------")
        print(f"Filtering Generation {pickaxe.generation}\n")

    def _pre_print(self) -> None:
        """
        Prints a message indicating that the filter is being applied.
        """
        print(f"Applying filter: {self.filter_name}")

    def _post_print(self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float) -> None:
        """
        Prints the results of the filtering process after it completes.

        Args:
            pickaxe (Pickaxe): The object representing the current pickaxe with reactions and compounds.
            n_total (int): The total number of compounds before filtering.
            n_filtered (int): The number of compounds remaining after filtering.
            time_sample (float): The time taken for the filtering process.
        """
        print(
            f"{n_filtered} of {n_total} compounds remain after applying "
            f"filter: {self.filter_name}"
            f"--took {round(time.time() - time_sample, 2)}s.\n"
        )

    def _post_print_footer(self, pickaxe: Pickaxe) -> None:
        """
        Prints a footer message after the filtering process is complete.

        Args:
            pickaxe (Pickaxe): The object representing the current pickaxe with reactions and compounds.
        """
        print(f"Done filtering Generation {pickaxe.generation}")
        print("----------------------------------------\n")

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int) -> Tuple[Set[str], Set[str]]:
        """
        Filters reactions based on their feasibility score and removes compounds 
        associated with unfeasible reactions.

        Args:
            pickaxe (Pickaxe): The object representing the current pickaxe with reactions and compounds.
            processes (int): The number of processes to use for parallelization (if applicable).
        
        Returns:
            Tuple[Set[str], Set[str]]: A tuple containing two sets:
                - 'reactions_to_remove': The reactions that have been filtered out.
                - 'compounds_to_remove' : The compounds associated with the removed reactions.
        """
        reactions_to_remove, compounds_to_remove = set(), set()

        if not pickaxe.reactions:
            print("No reactions found in pickaxe. Skipping filtering.")
            return reactions_to_remove, compounds_to_remove

        def get_cpd_smiles(cpd_id: str) -> str:
            """
            Helper function to retrieve the SMILES string of a compound by its ID.

            Args:
                cpd_id (str): The ID of the compound.

            Returns:
                str: The SMILES string of the compound.
            """
            return pickaxe.compounds.get(cpd_id, {}).get("SMILES", "")

        def _filter_reactions_by_feasibility() -> None:
            """
            Filters reactions based on their feasibility score. Removes reactions with scores 
            below the threshold and deletes compounds associated with those reactions.
            """
            current_generation = pickaxe.generation

            reactions_to_delete = set()
            compounds_to_remove = set()

            for rxn_id, reaction_data in list(pickaxe.reactions.items()):
                reactant_smiles = [get_cpd_smiles(cpd[1]) for cpd in reaction_data.get("Reactants", [])]
                product_smiles = [get_cpd_smiles(cpd[1]) for cpd in reaction_data.get("Products", [])]
                reaction_str = f"{' + '.join(reactant_smiles)} = {' + '.join(product_smiles)}"

                try:
                    feasibility_score = self.model.predict_proba(reaction_str)
                    feasibility_score = np.array(feasibility_score, dtype=float).flatten()
                except Exception:
                    feasibility_score = np.array([0.0])

                if np.any(feasibility_score < self.threshold):
                    reactions_to_delete.add(rxn_id)

                    for cpd in reaction_data.get("Reactants", []) + reaction_data.get("Products", []):
                        cpd_id = cpd[1]
                        if (
                            cpd_id in pickaxe.compounds
                            and pickaxe.compounds[cpd_id].get("Type") not in ["Coreactant", "Target Compounds", "Starting Compound"]
                        ):
                            compounds_to_remove.add(cpd_id)

            for rxn_id in reactions_to_delete:
                del pickaxe.reactions[rxn_id]

            for cpd_id in compounds_to_remove:
                del pickaxe.compounds[cpd_id]

            print(f"Deleted {len(reactions_to_delete)} reactions with feasibility scores below {self.threshold}.")
            print(f"Deleted {len(compounds_to_remove)} compounds associated with removed reactions.")

        _filter_reactions_by_feasibility()

        for cpd_id, cpd_data in pickaxe.compounds.items():
            if cpd_data.get("Type") in ["Coreactant", "Target Compounds", "Starting Compound"] or cpd_id not in compounds_to_remove:
                cpd_data["Expand"] = True

        return reactions_to_remove, compounds_to_remove

if __name__ == "__main__":
    pass
