'''
Define classes for pathway and reaction entries
'''

from src.utils import sort_x_by_y
from src.pathway_utils import get_stoich_pk
from collections import namedtuple, defaultdict
from typing import List

DatabaseEntry = namedtuple("DatabaseEntry", "db, id")
Enzyme = namedtuple("Enzyme", "uniprot_id, sequence")

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


class pathway:
    def __init__(self, rhashes, starter_hash=None, target_hash=None, prc_mcs=None):
        self.starter = starter_hash
        self.target = target_hash
        self.rhashes = rhashes # Hash ids for the path's reactions, in order
        self.prc_mcs = prc_mcs # Average (over known rxns) peri-rxn-ctr MCS score for each predicted rxn (tuple)
        self.mdf = None # Min-max driving force
        self.dG_opt = None # dGrs given optimized concentrations
        self.dG_err = None # Uncertainty about dG_opt
        self.conc_opt = None # Optimized substrate concentrations
        

    def min_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return min(self.prc_mcs)
        
    def max_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return max(self.prc_mcs)
        
    def mean_mcs(self):
        if self.prc_mcs is None:
            return None
        else:
            return sum(self.prc_mcs) / len(self.prc_mcs)

    def compute_mean_prc_mcs(self, pred_rxns):
        '''
        Assumes known reactions sorted by 
        substrate_averaged_prc_mcs and could
        break. Need to ensure entries always sorted
        in expansion object.
        '''
        # Computes an average of averages. First average
        # is over prc mcs of reactants in each known reaction
        # Second average over known reactions associated w/ 
        # a predicted reaction
        self.prc_mcs = []
        for rh in self.rhashes:
            krs = pred_rxns[rh].known_rxns
            kr_mean_mcs = 0
            for i, elt in enumerate(krs):
                if elt[0] is None:
                    break # This is where assumption of sorting comes in
                elif i > 0:
                    kr_mean_mcs = (kr_mean_mcs * i + sum(elt[0]) / len(elt[0])) / (i + 1) # Rolling ave
                else:
                    kr_mean_mcs = sum(elt[0]) / len(elt[0])

            self.prc_mcs.append(kr_mean_mcs)

class reaction:
    def __init__(self, rid, smarts, smi2pkid=None, rules=[], known_rxns=[]):
        self.rid = rid
        self.smarts = smarts
        self.rules = rules
        self.known_rxns = known_rxns
        self.smi2pkid = smi2pkid # pk ids : smiles

    def sort_known_rxns(self):
        '''
        Sort by mean for now. May include input to sort by min
        '''
        krs_w_mcs = [elt for elt in self.known_rxns if elt[0] is not None]
        krs_wo_mcs = [elt for elt in self.known_rxns if elt[0] is None]

        if krs_w_mcs:
            mcses = list(zip(*krs_w_mcs))[0]
            mean_mcses = list(map(lambda x: sum(x) / len(x), mcses))
            krs_w_mcs, _ = sort_x_by_y(krs_w_mcs, mean_mcses, reverse=True)
            self.known_rxns = list(krs_w_mcs) + krs_wo_mcs

# class expansion:
#     def __init__(self, paths=[], pred_rxns={}):
#         self.paths = paths
#         self.pred_rxns = pred_rxns

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

def get_smi2pkid(rhash, pk):
    smi2pkid = {}
    for elt in pk.reactions[rhash]['Reactants']:
        smi2pkid[pk.compounds[elt[1]]['SMILES']] = elt[1]

    for elt in pk.reactions[rhash]['Products']:
        smi2pkid[pk.compounds[elt[1]]['SMILES']] = elt[1]
    
    return smi2pkid