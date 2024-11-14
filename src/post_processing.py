from src.utils import load_json
from src.cheminfo_utils import standardize_smiles
from minedatabase.utils import get_compound_hash
from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict, Iterable, Any
from enum import Enum
import hashlib
import pathlib
import pickle
from rdkit import rdBase
from rdkit.Chem import CanonSmiles
from collections import defaultdict
from itertools import permutations, product
import networkx as nx
import re

# TODO: modify pickaxe so that we don't have to do this downstream
def realign_pred_rxn_to_rule(rxn_smarts: str, rule_template: str, coreactants: dict[str, str]) -> list[tuple[int]]:
    '''
    Returns permutations of reaction's reactant indices that match rule template

    Args
    ----
    rxn_smarts:str
        Reaction smarts
    rule_template:str
        Reactant names from rule ordered as designated in the operator.
        Separated by ';'
    coreactants:dict[str, str]
        Lookup SMILES to coreactant names

    Returns
    ------
    matched_idxs:list[tuple[int]]
        List of correctly ordered reactant indices
    '''
    rct_smi = rxn_smarts.split('>>')[0].split('.')
    rule_template = rule_template.split(';')

    if len(rct_smi) != len(rule_template):
        return []
    
    rct_names = [coreactants.get(smi, ('Any',)) for smi in rct_smi]
    rct_idxs = [i for i in range(len(rct_names))]
    matched_idxs = []
    for perm in permutations(rct_idxs):
        permuted_rct_names = [rct_names[idx] for idx in perm]
        for prod in product(*permuted_rct_names):
            if all([elt[0] == elt[1] for elt in zip(rule_template, prod)]):
                matched_idxs.append(perm)
        
    return matched_idxs

def get_path_id(reaction_ids:Iterable[str]):
    '''Returns hash id for pathway given
    reaction ids'''
    concat = "".join(reaction_ids)
    return "P" + hashlib.sha1(concat.encode("utf-8")).hexdigest()

@dataclass
class DatabaseEntry:
    name:str
    id:str

    @classmethod
    def from_dict(cls, dbe:dict):
        return cls(**dbe)

class Expansion:
    def __init__(self, forward: pathlib.Path = None, reverse: pathlib.Path = None, operator_reverses: dict = None):
        self.operator_reverses = operator_reverses

        if not forward and not reverse:
            raise ValueError("Must provide at least one expansion")
        
        if reverse and not operator_reverses:
            raise ValueError("Must provide mappings of operator to reverses to process a reverse expansion")
        
        self.forward = self._load(forward, flip=False) if forward else None
        self.reverse = self._load(reverse, flip=True) if reverse else None
        
        if self.forward and self.reverse:
            self.checkpoints = {k for k in self.forward['compounds'] if k[0] != 'X'} & {k for k in self.reverse['compounds'] if k[0] != 'X'}
        else:
            self.checkpoints = None

        self._compounds = None
        self._reactions = None
        self._coreactants = None

        if self.forward and self.reverse:
            self._compounds = {**self.forward['compounds'], **self.reverse['compounds']}
            self._reactions = {**self.forward['reactions'], **self.reverse['reactions']}
            self._coreactants = {**self.forward['coreactants'], **self.reverse['coreactants']}

    @property
    def coreactants(self):
        if self._coreactants:
            return self._coreactants
        elif not self.reverse:
            return self.forward['coreactants']
        else:
            return self.reverse['coreactants']
    
    @property
    def compounds(self):
        if self._compounds:
            return self._compounds
        elif not self.reverse:
            return self.forward['compounds']
        else:
            return self.reverse['compounds']
        
    @compounds.setter
    def compounds(self, value: dict):
        self._compounds = value
    
    @property
    def starters(self):
        if self.forward: # Forward or combo
            return self.forward['starters']
        else: # Retro
            return self.reverse['starters']
        
    @property
    def targets(self):
        if self.reverse: # Retro or combo
            return self.reverse['targets']
        else: # Forward
            return self.forward['targets']
        
    @property
    def reactions(self):
        if self._reactions:
            return self._reactions
        elif not self.reverse:
            return self.forward['reactions']
        else:
            return self.reverse['reactions']
        
    @reactions.setter
    def reactions(self, value: dict):
        self._reactions = value
        
    @property
    def generations(self):
        if self.forward and self.reverse:
            return self.forward['generations'] + self.reverse['generations']
        elif self.forward:
            return self.forward['generations']
        elif self.reverse:
            return self.reverse['generations']

    def _load(self, filepath: pathlib.Path, flip: bool) -> dict[str, Any]:
        '''
        Loads half expansions and stores in the direction of real-world synthesis

        Args
        -----
        filepath:pathlib.Path
            To half expansion file
        flip:bool
            Need to flip reactions, i.e., direction of expansion
            is opposite real-world direction of synthesis
        '''

        attributes = [
            'compounds',
            'reactions',
            'operators',
            'coreactants',
            'starters',
            'targets',
            'generations'
        ]
        half_expansion = {}

        with open(filepath, 'rb') as f:
            contents = pickle.load(f)

        for k in attributes:
            if k == 'generations':
                half_expansion[k] = contents[k]
            
            elif k == 'coreactants':
                coreactants = defaultdict(set)
                for name, (smiles, cid) in contents['coreactants'].items():
                    smiles = standardize_smiles(smiles, do_remove_stereo=True, do_find_parent=False, do_neutralize=False, quiet=True)
                    coreactants[smiles].add(name)
                    half_expansion[k] = coreactants
            
            elif k == 'starters':
                starters = {}
                for v in contents['compounds'].values():
                    if v["Type"].startswith("Start"):
                        starters[v['_id']] = v["ID"]

            elif k == 'targets':
                targets = {}
                for v in contents['targets'].values():
                    target_cid = get_compound_hash(v['SMILES'])[0] # Get neutral (non-target) cpd hash
                    target_name = v['ID']
                    targets[target_cid] = target_name

            elif flip and k == 'reactions':
                reversed_reactions = dict([self._flip_reaction(rxn, self.operator_reverses) for rxn in contents[k].values()])
                half_expansion[k] = reversed_reactions
            
            else:
                half_expansion[k] = contents[k]

        if flip:
            half_expansion['targets'] = starters
            half_expansion['starters'] = targets
        else:
            half_expansion['starters'] = starters
            half_expansion['targets'] = targets

        return half_expansion

    def _flip_reaction(self, rxn: dict, operator_reverses: dict):
        flipped_rxn = {}
        flipped_rxn['_id'] = get_reaction_hash(rxn['Products'], rxn['Reactants']) # TODO go with this when ready
        # flipped_rxn['_id'] = rxn['_id'] + '_reverse'
        flipped_rxn['Reactants'] = rxn['Products']
        flipped_rxn['Products'] = rxn['Reactants']
        flipped_operators = set()
        for o in rxn['Operators']:
            if o in operator_reverses:
                for fo in operator_reverses[o]:
                    flipped_operators.add(fo)
        flipped_rxn['Operators'] = flipped_operators
        flipped_rxn['SMILES_rxn'] = ' => '.join(rxn['SMILES_rxn'].split(' => ')[::-1])
        return flipped_rxn['_id'], flipped_rxn
    
    def find_paths(self):
        '''
        Returns reaction paths from real-world starters to real-world targets
        in real-world forward direction.
        '''
        def dictize_paths(paths: list):
            dpaths = defaultdict(list)
            for p in paths:
                dpaths[(p[0], p[-1])].append(p[1::2])
            return dpaths

        # Single direction expansions
        if self.forward and not self.reverse:
            paths = self._find_paths(self.forward, retro=False)
        elif self.reverse and not self.forward:
            paths = self._find_paths(self.reverse, retro=True)

        else: # Combo expansions
            forward_paths = self._find_paths(self.forward, retro=False)
            retro_paths = self._find_paths(self.reverse, retro=True)

            # Stitch together halves of combo
            forward_by_checkpoint = defaultdict(list)
            for fp in forward_paths:
                forward_by_checkpoint[fp[-1]].append(fp)

            retro_by_checkpoint = defaultdict(list)
            for rp in retro_paths:
                retro_by_checkpoint[rp[0]].append(rp)

            paths = []
            for ckpt in forward_by_checkpoint.keys() & retro_by_checkpoint.keys():
                path_pairs = product(forward_by_checkpoint[ckpt], retro_by_checkpoint[ckpt])
                for pair in path_pairs:
                    paths.append(pair[0][:-1] + pair[1])

            # Catch paths that connect over just half of the combo expansion
            for rws in self.forward['starters']:
                for p in retro_by_checkpoint[rws]:
                    paths.append(p)
            for rwt in self.reverse['targets']:
                for p in forward_by_checkpoint[rwt]:
                    paths.append(p)

        return dictize_paths(paths)     

    def _find_paths(self, half_expansion: dict, retro: bool):
        '''
        
        Args
        ----
        retro
            is retro expansion
        
        '''
        if self.checkpoints and retro: # Combo, retro half
            sources = self.checkpoints
            sinks = half_expansion['targets'].keys()
        elif self.checkpoints and not retro: # Combo, fwd half
            sources = half_expansion['starters'].keys()
            sinks = self.checkpoints
        else: # Single direction expansion, both stored in rw synth dir
            sources = half_expansion['starters'].keys()
            sinks = half_expansion['targets'].keys()

        # Faster to search from few to many
        if len(sources) > len(sinks):
            tmp = sources
            sources = sinks
            sinks = tmp
            flip = True
        else:
            flip = False

        sinks = list(sinks - sources) # all_simple_paths returns empty generator if source in sink

        DG = self.construct_network(half_expansion)

        if flip:
            DG = DG.reverse(copy=False)

        paths = []
        for sr in sources:
            paths += list(nx.all_simple_paths(DG, source=sr, target=sinks, cutoff=half_expansion['generations'] * 2)) # Search up to gens x 2 (bc bipartite)

        return [elt[::-1] for elt in paths] if flip else paths
    
    def construct_network(self, half_expansion: dict):
        '''
        Constructs bipartite (compounds and reactions) directed graph
        in the order of real-world synthesis

        Args
        ----
        half_expansion
        flip
            Flips orientation of the edges (direction of the reactions)
        '''
        cpd_node_list = []
        rxn_node_list = []
        edge_list = []

        # Get compound information
        for i, cpd in half_expansion['compounds'].items():

            if i[0] == "X": # Skip coreactants
                 continue

            cpd_node_list.append(
                (
                    i,
                    {
                        "SMILES": cpd["SMILES"],
                        "InChIKey": cpd["InChI_key"],
                        "Type": cpd["Type"],
                    }
                )
            )

        # Get reaction information
        for i, rxn in half_expansion['reactions'].items():

            rxn_node_list.append(
                (
                    i,
                    {
                        "Rule": rxn["Operators"],
                        "Type": "Reaction",
                    }
                )
            )

            # cpd_node_i => rxn_node_j only when cpd_node_i is only requirement beyond asssumed sources
            non_co_reactants = [c_id for _, c_id in rxn["Reactants"] if c_id[0] != 'X'] # Each entry in Reactants a unique molecule so list is ok
            if len(non_co_reactants) == 1:
                edge_list.append((non_co_reactants[0], i))

            for _, c_id in rxn["Products"]:

                if c_id[0] == 'X': # Don't bother w/ 'X' co-products
                     continue
                
                edge_list.append((i, c_id))
        
        # Create Graph
        DG = nx.DiGraph()  
        DG.add_nodes_from(cpd_node_list, bipartite=0)
        DG.add_nodes_from(rxn_node_list, bipartite=1)
        DG.add_edges_from(edge_list)
        
        return DG
    
    def prune(self, paths: dict[tuple, tuple]):
        '''
        Prune compounds and reactions to save time
        e.g., with thermo stuff
        '''
        pruned_rxns = {}
        pruned_cpds = {}

        for st_paths in paths.values():
            for p in st_paths:
                for rid in p:
                    rxn = self.reactions[rid]
                    pruned_rxns[rid] = rxn
                    for _, cid in rxn['Reactants'] + rxn['Products']:
                        pruned_cpds[cid] = self.compounds[cid]
        
        self.reactions = pruned_rxns
        self.compounds = pruned_cpds
    
def get_stoich_pk(smiles_rxn):
    lhs, rhs = [r for r in smiles_rxn.split("=>")]
    lhs, rhs = [[rct.strip(" ") for rct in side.split(" + ")] for side in [lhs, rhs]]
    [rct.split(" ") for rct in lhs]
    pat = re.compile("\((\d+)\) (.+)")
    reactants = {get_canon_smiles(smiles): -1*int(stoich) for stoich, smiles in [pat.findall(rct)[0] for rct in lhs]}
    products = {get_canon_smiles(smiles): int(stoich) for stoich, smiles in [pat.findall(rct)[0] for rct in rhs]}
    reactants.update(products)
    return reactants

def get_canon_smiles(smi):
    _ = rdBase.BlockLogs()
    try:
        return CanonSmiles(smi)
    except BaseException as e:
        return smi
    
def get_reaction_hash(
    reactants: List[tuple], products: List[tuple]
) -> str:
    """Tries to do what Pickaxe does mid transform.
    By the way pickaxe is not doing what I think it
    wants to do.
    TODO: make hashes stoich-sensitive
    """
    def to_str(half_rxn):
        return [f"(1) {x}" for x in sorted(half_rxn)]

    reactant_ids = [reactant[1] for reactant in reactants]
    product_ids = [product[1] for product in products]
    text_ids_rxn = (
        " + ".join(to_str(reactant_ids)) + " => " + " + ".join(to_str(product_ids))
    )
    rhash = "R" + hashlib.sha256(text_ids_rxn.encode()).hexdigest()

    return rhash

@dataclass
class Enzyme:
    uniprot_id:str
    sequence:Optional[str] = None
    existence:Optional[str] = None
    reviewed:Optional[str] = None
    ec:Optional[str] = None
    organism:Optional[str] = None
    name:Optional[str] = None

    def to_dict(self):
        return asdict(self)
    
    @classmethod
    def from_dict(cls, enz:dict):
        return cls(**enz)

@dataclass
class KnownReaction:
    id:str
    smarts:str
    operators:List[str]
    enzymes:List[Enzyme]
    db_entries:List[DatabaseEntry]
    reaction_center:Iterable[Iterable] = field(default_factory=tuple)
    image:str = ''

    def to_dict(self):
        return asdict(self)
    
    @classmethod
    def from_dict(cls, kr:dict):
        '''
        Construct from dict

        Args
        ----
        kr:dict
            cls.to_dict() output
        '''
        objectifiers = {
            'enzymes': lambda L : [Enzyme.from_dict(e) for e in L],
            'db_entries': lambda L : [DatabaseEntry.from_dict(d) for d in L]
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in kr.items())
        return cls(**kwargs)
    
    def filter_enzymes(self, filter_by:Dict[str, Iterable]):
        '''
        Removes enzymes without desired value in given field

        Args
        ----
        filter_by:Dict[str, Iterable]
            Maps field name -> list of valid values
        '''
        
        tmp = []
        for e in self.enzymes:
            include = True
            for k, v in filter_by.items():
                if getattr(e, k) not in v:
                    include = False
                    break
            
            if include:
                tmp.append(e)
        
        self.enzymes = tmp 

    @property
    def has_valid_enzymes(self):
        return True if len(self.enzymes) > 0 else False

@dataclass
class PredictedReaction:
    id:str
    smarts:str
    operators:List[str]
    reaction_center:Iterable[Iterable] = field(default_factory=tuple)
    analogues:Dict[str, KnownReaction] = field(default_factory=dict)
    rcmcs:Dict[str, float] = field(default_factory=dict)
    image:str = ''

    def to_dict(self):
        '''Returns dict representation w/ all fields
        except for analogues, only ids are kept'''
        self_dict = asdict(self)
        self_dict['analogues'] = list(self_dict['analogues'].keys())
        return self_dict

    @classmethod
    def from_dict(cls, pr:dict, stored_known_reactions:Dict[str, KnownReaction]):
        '''
        Construct from stored predicted reaction and stored known reactions

        Args
        -----
        pr:dict
            cls.to_dict() output
        stored_known_reactions:dict
            Contains KnownReaction.to_dict() outputs stored under known 
            reaction ids as keys
        '''
        objectifiers = {
            'analogues': lambda L : {krid : stored_known_reactions[krid] for krid in L} if L else {}
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in pr.items())
        return cls(**kwargs)
    
    @classmethod
    def from_pickaxe(cls, pk, id):
        smarts = pk_rhash_to_smarts(id, pk)
        operators = list(pk.reactions[id]["Operators"])
        return cls(id, smarts, operators)
    
    def filter_analogues(self, filter_by:Dict[str, Iterable]):
        tmp = {}
        for krid, a in self.analogues.items():
            include = True
            for k, v in filter_by.items():
                if getattr(a, k) not in v:
                    include = False
                    break
            
            if include:
                tmp[krid] = a
        
        self.analogues = tmp
    
    @property
    def has_valid_analogues(self):
        if len(self.analogues) == 0:
            return False
        elif not all([a.has_valid_enzymes for a in self.analogues.values()]):
            return False
        else:
            return True
    
    def top_rcmcs(self, k:int):
        '''Top k RCMCS scores'''
        srt_krids = sorted(self.rcmcs, key=lambda x : self.rcmcs[x], reverse=True)[:k]
        return [self.rcmcs[krid] for krid in srt_krids]
    
    def top_analogues(self, k:int):
        '''Top k analogues, scored on RCMCS'''
        srt_krids = sorted(self.rcmcs, key=lambda x : self.rcmcs[x], reverse=True)[:k]
        return [self.analogues[krid] for krid in srt_krids]

@dataclass
class Path:
    id:str
    starter:str
    target:str
    reactions:List[PredictedReaction]
    mdf:float
    dG_opt:Dict[str, float]
    dG_err:Dict[str, float]
    sid:str # Starter hash
    tid:str # Target hash

    def to_dict(self):
        '''Returns dict representation w/ all fields
        except for reactions, only ids are kept'''
        self_dict = asdict(self)
        self_dict['reactions'] = [r.id for r in self.reactions]
        return self_dict

    @classmethod
    def from_dict(cls, path:dict, prs:Dict[str, PredictedReaction]):
        '''
        Construct from stored path, stored predicted reactions, &
        stored known reactions

        Args
        -----
        path:dict
            cls.to_dict() output
        stored_known_reactions:dict
            Contains PredictedReaction.to_dict() outputs stored under predicted 
            reaction ids as keys
        stored_known_reactions:dict
            Contains KnownReaction.to_dict() outputs stored under known 
            reaction ids as keys
        '''
        objectifiers = {
            'reactions': lambda L : [prs[prid] for prid in L] if L else []
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in path.items())
        return cls(**kwargs)
    
    @property
    def valid(self):
        return all([r.has_valid_analogues for r in self.reactions])
    
    @property
    def min_rcmcs(self):
        return min([r.top_rcmcs(k=1)[0] for r in self.reactions])
    
    @property
    def mean_rcmcs(self):
        top_rcmcs = [r.top_rcmcs(k=1)[0] for r in self.reactions]
        return sum(top_rcmcs) / len(top_rcmcs)
    
    def aggregate_rcmcs(self, pr_agg:str, kr_agg:str, k:int):
        aggs = {
            'min': lambda x : min(x),
            'mean': lambda x: sum(x) / len(x),
            'max': lambda x: max(x)
        }
        if pr_agg not in aggs or kr_agg not in aggs:
            raise ValueError(f"Choose valid aggregation methods from: {aggs.keys()}")
        
        top_k_rcmcs = [r.top_rcmcs(k=k) for r in self.reactions]
        kr_agged_rcmcs = [aggs[kr_agg](top) for top in top_k_rcmcs]
        return aggs[pr_agg](kr_agged_rcmcs)
    
class EnzymeExistence(Enum):
    PROTEIN = 'Evidence at protein level'
    TRANSCRIPT = 'Evidence at transcript level'
    HOMOLOGY = 'Inferred from homology'
    PREDICTED = 'Predicted'
    UNCERTAIN = 'Uncertain'

class PathWrangler:
    '''
    Loads dict representations of stored paths, 
    predicted reactions, & known reactions and provides
    method to get valid paths with options to filter & sort
    based on various criteria. 
    '''
    enzyme_existence = EnzymeExistence
    
    def __init__(self, proc_exp: pathlib.Path, img_subdir: str) -> None:
        img_prefix = proc_exp / img_subdir
        self.known_reactions = self._prepend_images(load_json(proc_exp / "known_reactions.json"), img_prefix)
        self.predicted_reactions = self._prepend_images(load_json(proc_exp / "predicted_reactions.json"), img_prefix)
        self.paths = load_json(proc_exp / "found_paths.json")
        self.starter_targets = self._extract_starter_targets()

    def _prepend_images(self, d: dict, prefix: pathlib.Path):
        for v in d.values():
            v['image'] = prefix / v['image']

        return d
        
    def _extract_starter_targets(self):
        tmp = set()
        for k, v in self.paths.items():
            tmp.add((v['starter'], v['target']))
        return tuple(tmp)

    def get_paths(
            self,
            starters:Iterable[str] = None,
            targets:Iterable[str] = None,
            filter_by_enzymes:Dict[str, Iterable] = None,
            sort_by= None,
        ) -> list[Path]:
        '''
        Get valid paths with options to filter & sort

        Args
        -------
        starters:Iterable[str]
            Names of desired starter molecules
        targets:Iterable[str]
            Names of desired target molecules
        filter_by_enzymes:Dict[str, Iterable]
            Maps field name -> list of valid values
        sort_by:str | dict
            'min_rcmcs' sorts by path's minimum reaction center mcs
                given by taking max RCMCS for each predicted reaction 
            'mean_rcmcs' sorts by paths' mean reaction center mcs
                given by taking max RCMCS for each predicted reaction
            You may provide a dict that specifies how to aggregate RCMCS on
                both the path level (multiple predicted reactions) and the 
                predicted reaction level (multiple analogues)
                In {'kr_agg': str, 'pr_agg': str, 'k': int}
                pr_agg aggregates on path level 'min' | 'mean' | 'max'
                kr_agg aggregates on predicted reaction level 'min' | 'mean' | 'max'
                k takes k top analogues per predicted reaction

        Returns
        ---------
        Fitlered & sorted List[Path]
        '''
        if starters and targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['starter'] in starters and v['target'] in targets]
        elif starters and not targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['starter'] in starters]
        elif not starters and targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['target'] in targets]
        else:
            req_path_ids = [v['id'] for v in self.paths.values()]

        if filter_by_enzymes:
            if 'existence' in filter_by_enzymes:
                filter_by_enzymes['existence'] = [self.enzyme_existence[v.upper()].value for v in filter_by_enzymes['existence']]
        
        req_prids = set()
        req_krids = set()
        for pid in req_path_ids:
            for prid in self.paths[pid]['reactions']:
                req_prids.add(prid)
                for krid in self.predicted_reactions[prid]['analogues']:
                    req_krids.add(krid)

        krs = {}
        for krid in req_krids:
            kr = KnownReaction.from_dict(self.known_reactions[krid])
            if filter_by_enzymes:
                kr.filter_enzymes(filter_by=filter_by_enzymes)
            krs[krid] = kr

        prs = {}
        for prid in req_prids:
            pr = PredictedReaction.from_dict(self.predicted_reactions[prid], krs)
            prs[prid] = pr

        paths = []
        for pid in req_path_ids:
            p = Path.from_dict(self.paths[pid], prs)
            if p.valid:
                paths.append(p)

        if sort_by is None:
            return paths
        elif type(sort_by) is str:
            return sorted(paths, key= lambda p : getattr(p, sort_by), reverse=True)
        elif type(sort_by) is dict:
            self._validate_custom_agg(sort_by)
            return sorted(paths, key= lambda p : p.aggregate_rcmcs(**sort_by))
        
    def get_path_with_id(self, pid: str) -> Path:
        '''
        Get path with id pid
        '''
        # Required reaction ids for this path id
        req_prids = set()
        req_krids = set()
        for prid in self.paths[pid]['reactions']:
            req_prids.add(prid)
            for krid in self.predicted_reactions[prid]['analogues']:
                req_krids.add(krid)

        # Get known reactions
        krs = {}
        for krid in req_krids:
            kr = KnownReaction.from_dict(self.known_reactions[krid])
            krs[krid] = kr

        # Get predicted reactions
        prs = {}
        for prid in req_prids:
            pr = PredictedReaction.from_dict(self.predicted_reactions[prid], krs)
            prs[prid] = pr

        p = Path.from_dict(self.paths[pid], prs) # Construct path

        if not p.valid:
            print("Warning, this path may not have a complete set of known analogues & enzymes!")

        return p
    
    def _validate_custom_agg(self, sort_by:dict):
        req_kwargs = {'kr_agg': str, 'pr_agg': str, 'k': int}
        for k, t in req_kwargs.items():
            if k not in sort_by:
                raise ValueError(f"Missing required argument {k}")
            if type(sort_by[k]) is not t:
                raise ValueError(f"Invalid type {sort_by[k]} for argument {k}. Must be {t}")
            
def pk_rhash_to_smarts(rhash: str, pk: Expansion):
    '''
    Make reaction smarts string for
    reaction indexed by rhash in a Expansion
    object
    '''
    rxn_smiles = pk.reactions[rhash]["SMILES_rxn"]
    rxn_stoich = get_stoich_pk(rxn_smiles)
    products = ".".join([".".join([smi]*stoich) for smi, stoich in rxn_stoich.items() if stoich >= 0])
    reactants = ".".join([".".join([smi]*abs(stoich)) for smi, stoich in rxn_stoich.items() if stoich <= 0])
    rxn_sma = ">>".join([reactants, products])
    return rxn_sma

if __name__ == '__main__':
    # e1 = Enzyme('up1', 'AARTGW', 'Evidence at protein level', 'reviewed', '1.1.1.1', 'mouse', 'mouse protein')
    # e2 = Enzyme('up2', 'QWYPOI', 'Evidence at transcript level', 'reviewed', '1.1.1.2', 'e. coli', 'bacteria protein')
    # db1 = DatabaseEntry('rhea', '123')
    # db2 = DatabaseEntry('rhea', '456')
    # es = [e1, e2]
    # dbs = [db1, db2]
    # kr = KnownReaction('1', 'CC=O>>CCO', '/home/stef/bottle/artifacts/imgs/img1.svg', es, dbs)
    # kr2 = KnownReaction('2', 'CC(=O)O>>CC.O=C=O', '/home/stef/bottle/artifacts/imgs/img2.svg', es, dbs)
    # pr = PredictedReaction(
    #     id='1',
    #     smarts='O=CCC(=O)O>>O=CCC.O=C=O',
    #     image='/home/stef/bottle/artifacts/imgs/img3.svg',
    #     operators=['rule1', 'rule2'],
    #     analogues={'1':kr, '2': kr2},
    #     rcmcs={'1': 0.9, '2': 0.2}
    # )

    # # to / from dict
    # kr_dict = asdict(kr)
    # kr_from = kr.from_dict(kr_dict)
    # pr_dict = pr.to_dict()
    # pr_from = PredictedReaction.from_dict(pr_dict, {'1':kr, '2':kr2})

    # # Filter enzymes in known reactions
    # kr.filter_enzymes({'existence': [EnzymeExistence.PROTEIN.value]})
    # print(kr.has_valid_enzymes)
    # kr2.filter_enzymes({'existence': [EnzymeExistence.HOMOLOGY.value]})

    # # Filter known reactions in predicted reaction
    # pr.filter_analogues({'has_valid_enzymes': [True]})

    # # Sort analogues
    # pr = PredictedReaction(
    #     id='1',
    #     smarts='O=CCC(=O)O>>O=CCC.O=C=O',
    #     image='/home/stef/bottle/artifacts/imgs/img3.svg',
    #     operators=['rule1', 'rule2'],
    #     analogues={'1':kr, '2': kr2},
    #     rcmcs={'1': 0.9, '2': 0.2}
    # )

    # print(pr.top_analogues(k=2))
    # print(pr.top_analogues(k=1))

    # # Top RCMCS
    # print(pr.top_rcmcs(k=2))
    # print(pr.top_rcmcs(k=1))

    # # Path-level operations
    # path = Path(
    #     id='1',
    #     starter='akg',
    #     target='hopa',
    #     reactions=[pr, pr],
    #     sid='Csdlfk',
    #     tid="Caweor34np"
    # )

    # print(path.valid)
    # print(path.min_rcmcs)
    # print(path.mean_rcmcs)
    # print(path.aggregate_rcmcs(pr_agg='min', kr_agg='mean', k=5))

    # # Path wrangling
    # path_filepath = '../artifacts/processed_expansions/found_paths.json'
    # predicted_reactions_filepath = "../artifacts/processed_expansions/predicted_reactions.json"
    # known_reactions_filepath = "../artifacts/processed_expansions/known_reactions.json"
    # pw = PathWrangler(
    #     path_filepath=path_filepath,
    #     pr_filepath=predicted_reactions_filepath,
    #     kr_filepath=known_reactions_filepath
    # )

    # pw.get_paths()
    # pw.get_paths(sort_by='min_rcmcs')
    # pw.get_paths(sort_by='mean_rcmcs')
    # pw.get_paths(filter_by_enzymes={'existence':['protein']})
    
    # print('hold')

    from src.config import filepaths
    forward = filepaths['raw_expansions'] / "2_steps_pivalic_acid_to_bottle_targets_24_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_True.pk"
    expansion = Expansion(forward=forward)
    paths = expansion.find_paths()

    forward = filepaths['raw_expansions'] / "1_steps_alpha_ketoglutarate_to_None_rules_JN3604IMT_rules_carbonyl_free_co_metacyc_coreactants_carbonyl_free_sampled_False_pruned_False.pk"
    reverse = filepaths['raw_expansions'] / "1_steps_hopa_to_None_rules_JN3604IMT_rules_carbonyl_free_co_metacyc_coreactants_carbonyl_free_sampled_False_pruned_False.pk"
    imt_reverses = load_json(filepaths['rules'] / "jnimt_reverses.json")
    expansion = Expansion(forward=forward, reverse=reverse, operator_reverses=imt_reverses)
    paths = expansion.find_paths()
    print()
