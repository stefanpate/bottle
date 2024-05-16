'''
Define classes for pathway and reaction entries
'''

from src.utils import sort_x_by_y
from src.pathway_utils import get_stoich_pk
from collections import namedtuple, defaultdict
from typing import List
import numpy as np
import json

DatabaseEntry = namedtuple("DatabaseEntry", "db, id", defaults=[None, None])
Enzyme = namedtuple("Enzyme", "uniprot_id, sequence, ec, validation_score, existence, reviewed, organism", defaults=[None, None, None, None, None, None, None])

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
        self._sort_enzymes()
        self.enzyme_max_val = self.enzymes[0].validation_score if len(enzymes) > 0 else 0

    def _sort_enzymes(self):
        self.enzymes.sort(key=lambda e : e.validation_score, reverse=True)

class PredictedReaction(Reaction):
    _sort_keys = {'prc_mcs_min': lambda mcs_analogue: min(mcs_analogue[0]),
                    'prc_mcs_mean': lambda mcs_analogue: sum(mcs_analogue[0]) / len(mcs_analogue[0]),
                    'enzyme_validation': lambda mcs_analogue: mcs_analogue[1].enzyme_max_val}

    def __init__(self, id, smarts, imt_rules=[], analogues:List[KnownReaction]=[]):
        super().__init__(id, smarts, imt_rules)
        # TODO: ensure no repeat analogues
        prc_mcs_vals = [[0] for _ in range(len(analogues))]
        self._mcs_analogues = [list(elt) for elt in zip(prc_mcs_vals, analogues)]
        self.analogues_sorted_by = None
        self.substrates_reduced_by = None
        self.analogues = [elt[1] for elt in self._mcs_analogues]
        self._max_vals = {'prc_mcs_min': None, 'prc_mcs_mean': None, 'enzyme_validation':None}

    def max_prc_mcs(self, reduce_substrates):
        max_key = f"prc_mcs_{reduce_substrates}"
        if self._max_vals[max_key] is not None:
            return self._max_vals[max_key]
        else:
            self.sort_analogues([max_key], reduce_substrates)
            if reduce_substrates == 'min':
                max_val = min(self._mcs_analogues[0][0])
            elif reduce_substrates == 'mean':
                max_val = sum(self._mcs_analogues[0][0]) / len(self._mcs_analogues[0][0])
            self._max_vals[max_key] = max_val
            
            return max_val
        
    def max_enzyme_validation(self):
        max_key = "enzyme_validation"
        if self._max_vals[max_key] is not None:
            return self._max_vals[max_key]
        elif not self._mcs_analogues: # No analogues for this pr
            return 0
        else:
            self.sort_analogues([max_key], reduce_substrates=None)
            max_val = self._mcs_analogues[0][1].enzyme_max_val
            self._max_vals[max_key] = max_val
            return max_val

    def _get_sort_key(self, sort_by, reduce_substrates):
        fs = []
        for elt in sort_by:
            if elt == 'prc_mcs':
                fs.append(PredictedReaction._sort_keys[f"{elt}_{reduce_substrates}"])
            else:
                fs.append(PredictedReaction._sort_keys[elt])

        return lambda mcs_analogue : tuple(f(mcs_analogue) for f in fs)

    def sort_analogues(self, sort_by, reduce_substrates):
        if self.analogues_sorted_by == sort_by and self.substrates_reduced_by == reduce_substrates:
            pass
        else:
            sort_key = self._get_sort_key(sort_by, reduce_substrates)
            self._mcs_analogues.sort(reverse=True, key=sort_key)
            self.analogues = [elt[1] for elt in self._mcs_analogues]
            self.analogues_sorted_by = sort_by
            self.substrates_reduced_by = reduce_substrates

    def threshold(self, filter_by, reduce_substrates):
        bool_list = []
        for criterion, val in filter_by.items():
            if criterion == 'prc_mcs':
                max_val = self.max_prc_mcs(reduce_substrates)
            elif criterion == 'enzyme_validation':
                max_val = self.max_enzyme_validation()

            bool_list.append(max_val >= val)

        return all(bool_list)
    
    def top_analogue(self):
        if self._mcs_analogues:
            prc_mcs, analogue = self._mcs_analogues[0]
            
            # Reduce by whatever way was sorted by
            # Include indication of what way that was
            # in res dict
            if self.substrates_reduced_by == 'mean':
                prc_mcs = sum(prc_mcs) / len(prc_mcs)
            elif self.substrates_reduced_by == 'min':
                prc_mcs = min(prc_mcs)

            res = {'analogue':analogue, 'prc_mcs':prc_mcs, 'enzyme_validation':analogue.enzyme_max_val,
                    'sort_by':self.analogues_sorted_by, 'reduce_substrates':self.substrates_reduced_by}
            return res
        else:
            return {'analogue':None, 'prc_mcs':0, 'enzyme_validation':0,
                    'sort_by':self.analogues_sorted_by, 'reduce_substrates':self.substrates_reduced_by}
    
class Path:
    def __init__(self, starter, target, reaction_ids:List[str], id=None):
        self.id = id
        self.reaction_ids = reaction_ids
        self.starter = starter
        self.target = target
        self.prc_mcs_vals = [] # Average (over known rxns) peri-rxn-ctr MCS score for each predicted rxn (tuple)
        self.mdf = -np.inf # Min-max driving force
        self.dG_opt = [] # dGrs given optimized concentrations
        self.dG_err = [] # Uncertainty about dG_opt
        self.conc_opt = [] # Optimized substrate concentrations


class ProcessedExpansion:
    def __init__(self) -> None:
        self.starter_target_pairs = set()
        self._st2paths = defaultdict(list)
        self._id2path = {}
        self.predicted_reactions = {}
        self._next_path_id = 1

    def get_path_prc_mcs(self, path):
        # TODO: make it more obvious how prs are sorted
        # maybe make top_analogues take reduce substrate parameter optional
        # if None, uses stored reducing criteria (and sorting criteria)
        # and returns this...
        prc_mcs = []
        for rid in path.reaction_ids:
            pr = self.predicted_reactions[rid]
            prc_mcs.append(pr.top_analogue()['prc_mcs'])

        return prc_mcs
    
    def get_path_enzymes(self, path):
        enzymes = []
        for rid in path.reaction_ids:
            pr = self.predicted_reactions[rid]
            enzymes.append(pr.top_analogue()['analogue'].enzymes)

        return enzymes

    def add_path(self, starter, target, predicted_reactions:List[PredictedReaction]):
        # Add path to st2paths dict
        reaction_ids = [elt.id for elt in predicted_reactions]
        path = Path(starter, target, reaction_ids, self._next_path_id)
        self._st2paths[(starter, target)].append(path)
        self._id2path[self._next_path_id] = path
        self._next_path_id += 1

        # Add all new predicted reactions
        for pred_rxn in predicted_reactions:
            rid = pred_rxn.id
            if rid not in self.predicted_reactions:
                self.predicted_reactions[rid] = pred_rxn

        # Add new starter target pairs
        self.starter_target_pairs.add((starter, target))

    def get_paths_w_id(self, ids, sort_by=[], filter_by={}, reduce_substrates='mean', reduce_predicted_reactions='mean'):
        if type(ids) is int:
            ids = [ids]

        paths = [self._id2path[id] for id in ids]
        if not sort_by and not filter_by:
            return paths
        else:
            return self._get_paths(paths, sort_by, filter_by, reduce_substrates, reduce_predicted_reactions)
        
    def get_paths_w_st(self, starter:str, target:str, sort_by=[], filter_by={}, reduce_substrates='mean', reduce_predicted_reactions='mean'):
        paths = self._st2paths[(starter, target)] # Get paths for this starter-target pair
        if not sort_by and not filter_by:
            return paths
        else:
            return self._get_paths(paths, sort_by, filter_by, reduce_substrates, reduce_predicted_reactions)

    def _get_paths(self, paths, sort_by, filter_by, reduce_substrates, reduce_predicted_reactions):
        path_level_filters = ('mdf', )
        rxn_level_filters = ('enzyme_validation', 'prc_mcs')
        path_level_filter_by = {k : v for k, v in filter_by.items() if k in path_level_filters}
        rxn_level_filter_by = {k : v for k, v in filter_by.items() if k in rxn_level_filters}

        # Filter paths by path-level criteria
        if path_level_filter_by:
            for criterion, threshold in path_level_filter_by.items():
                paths = set(filter(lambda path : getattr(path, criterion) >= threshold, paths))
        
        # Get all unique predicted reactions for this selection
        prs = set()
        for path in paths:
            for rid in path.reaction_ids:
                prs.add(rid)
        
        # Filter by rxn-level criteria
        if rxn_level_filter_by:
            # Check which predicted reactions above thresholds
            pr_above_threshold = {}
            for rid in prs:
                pr_above_threshold[rid] = self.predicted_reactions[rid].threshold(rxn_level_filter_by, reduce_substrates)
            
            paths = set(filter(lambda path : all([pr_above_threshold[rid] for rid in path.reaction_ids]), paths)) # Filter paths

            prs = [rid for rid in prs if pr_above_threshold[rid]] # Filter pred reactions
        
        else:
            prs = list(prs)

        # Calculate measures for prs required by sort & filter
        pr_to_sort_key = {}
        for rid in prs:
            self.predicted_reactions[rid].sort_analogues(sort_by, reduce_substrates)
            top_analogue = self.predicted_reactions[rid].top_analogue()
            sort_key = [top_analogue[criterion] for criterion in sort_by]
            pr_to_sort_key[rid] = sort_key

        # Reduce sort keys along pred reactions in path
        path_id_to_sort_key = {}
        for path in paths:
            split_sort_keys = list(zip(*[pr_to_sort_key[rid] for rid in path.reaction_ids]))

            if reduce_predicted_reactions == 'min':
                path_id_to_sort_key[path.id] = tuple([min(elt) for elt in split_sort_keys])
            elif reduce_predicted_reactions == 'mean':
                path_id_to_sort_key[path.id] = tuple([sum(elt) / len(elt) for elt in split_sort_keys])

        # Sort paths
        sorted_ids = sorted(path_id_to_sort_key.keys(), reverse=True, key=lambda pid : path_id_to_sort_key[pid])
        sorted_paths = self.get_paths_w_id(sorted_ids)
        return sorted_paths
    
def load_known_rxns(path):
    with open(path, 'r') as f:
        data = json.load(f)

    for _,v in data.items():

        # Convert enzymes and db entries to namedtuples
        enzymes = []
        for e in v['enzymes']:
            for i in range(len(e)):
                if type(e[i]) == list: # Convert list ec to tuple for hashing, set ops
                    e[i]= tuple(e[i])

            enzymes.append(Enzyme(*e))


        # enzymes = [Enzyme(*elt) for elt in v['enzymes']]

        # # Convert EC number list to tuple
        # for e in enzymes:
        #     e.ec = tuple(e.ec)

        v['enzymes'] = enzymes

        db_entries = [DatabaseEntry(*elt) for elt in v['db_entries']]
        v['db_entries'] = db_entries

    return data

if __name__ == '__main__':
    import pickle

    # Params
    starters = 'ccm_v0'
    targets = 'hopa'
    generations = 3

    expansion_dir = '../data/processed_expansions/'
    fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pkl" # Expansion file name

    # Load processed expansions
    with open(expansion_dir + fn, 'rb') as f:
        pe = pickle.load(f)

    to_vis = {}
    for st_pair in pe.starter_target_pairs:
        starter, target = st_pair
        enzyme_validation_threshold = 0.95
        sort_by = ['enzyme_validation', 'prc_mcs']
        filter_by = {'mdf':0, 'enzyme_validation':enzyme_validation_threshold}

        paths = pe.get_paths_w_st(starter=starter,
                        target=target,
                        sort_by=sort_by,
                        filter_by=filter_by,
                        reduce_predicted_reactions='min'
                        )
    
    print('ok')