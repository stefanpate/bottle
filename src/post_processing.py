'''
Define classes for pathway and reaction entries
'''

from src.utils import sort_x_by_y
from src.pathway_utils import get_stoich_pk
from collections import namedtuple, defaultdict
from typing import List

DatabaseEntry = namedtuple("DatabaseEntry", "db, id", defaults=[None, None])
Enzyme = namedtuple("Enzyme", "uniprot_id, sequence, validation_score, existence, reviewed, organism", defaults=[None, None, None, None, None, None])

class Reaction:
        def __init__(self, id, smarts, imt_rules=[]):
            self.id = id
            self.smarts = smarts
            self.imt_rules = imt_rules

class KnownReaction(Reaction):
    def __init__(self, id, smarts, imt_rules=[], database_entries:List[DatabaseEntry]=[], enzymes:List[Enzyme]=[]):
        super().__init__(id, smarts, imt_rules)
        self.database_entries = database_entries
        self.enzymes = enzymes


    # def ave_prc_mcs(self):
    #     return sum(self.prc_mcs_vals) / len(self.prc_mcs_vals)

class PredictedReaction(Reaction):
    def __init__(self, id, smarts, imt_rules=[], analogues:List[KnownReaction]=[]):
        super().__init__(id, smarts, imt_rules)
        # TODO: ensure no repeat analogues
        prc_mcs_vals = [None for i in range(len(analogues))]
        self._mcs_analogues = list(zip(prc_mcs_vals, analogues))
        self._mcs_analogues.sort(reverse=True, key=self._mcs_analogue_sort_key)
        self.analogues = [elt[1] for elt in self._mcs_analogues]

    def _mcs_analogue_sort_key(self, elt):
        if elt[0] is None:
            return 0
        else:
            return elt[0]

    # def add_analogues(self, analogues:list):
    #     self.analogues += analogues
    #     self.analogues = list(set(self.analogues))
    #     self._sort_analogues()

    # def _sort_analogues(self):
    #     self.analogues.sort(reverse=True, key=lambda x : x.ave_prc_mcs)

class Path:
    def __init__(self, starter, target, reaction_ids:List[str]):
        self.reaction_ids = reaction_ids
        self.starter = starter
        self.target = target
        self.prc_mcs_vals = [] # Average (over known rxns) peri-rxn-ctr MCS score for each predicted rxn (tuple)
        self.mdf = None # Min-max driving force
        self.dG_opt = [] # dGrs given optimized concentrations
        self.dG_err = [] # Uncertainty about dG_opt
        self.conc_opt = [] # Optimized substrate concentrations


class ProcessedExpansion:
    def __init__(self) -> None:
        self._st2paths = defaultdict(list)
        self.predicted_reactions = {}

    def add_path(self, starter, target, predicted_reactions:List[PredictedReaction]):
        reaction_ids = [elt.id for elt in predicted_reactions]
        path = Path(starter, target, reaction_ids)
        self._st2paths[(path.starter, path.target)].append(path)

        for pred_rxn in predicted_reactions:
            rid = pred_rxn.id
            if rid not in self.predicted_reactions:
                self.predicted_reactions[rid] = pred_rxn

        # TODO: Check for and remove duplicate paths in st2paths
        # TODO: Keep lists in st2paths ordered

def rxn_hash_2_rxn_sma(rhash, pk):
    '''
    Make reaction smarts string for
    reaction indexed by rhash in a pk
    object
    '''
    rxn_stoich = get_stoich_pk(rhash, pk)
    products = ".".join([".".join([smi]*stoich) for smi, stoich in rxn_stoich.items() if stoich >= 0])
    reactants = ".".join([".".join([smi]*abs(stoich)) for smi, stoich in rxn_stoich.items() if stoich <= 0])
    rxn_sma = ">>".join([reactants, products])
    return rxn_sma