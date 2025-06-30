from src.utils import load_json
from src.cheminfo_utils import standardize_smiles
from minedatabase.utils import get_compound_hash
from typing import Iterable, Any
from enum import Enum
import hashlib
import pathlib
import pickle
from rdkit import rdBase
from rdkit.Chem import CanonSmiles
from collections import defaultdict
from itertools import product
import networkx as nx
import re
import polars as pl

def get_path_id(reaction_ids:Iterable[str]):
    '''Returns hash id for pathway given
    reaction ids'''
    concat = "".join(reaction_ids)
    return "P" + hashlib.sha1(concat.encode("utf-8")).hexdigest()

class Expansion:
    def __init__(
        self,
        forward: pathlib.Path = None,
        reverse: pathlib.Path = None,
        override_starters: dict[str, str] = None,
        override_targets: dict[str, str] = None,
    ):

        if not forward and not reverse:
            raise ValueError("Must provide at least one expansion")
        
        self.forward = self._load(forward, flip=False) if forward else None
        self.reverse = self._load(reverse, flip=True) if reverse else None
        
        if self.forward and self.reverse:
            self.checkpoints = {k for k in self.forward['compounds'] if k[0] != 'X'} & {k for k in self.reverse['compounds'] if k[0] != 'X'}
            self.forward['targets'] = self.checkpoints
            self.reverse['starters'] = self.checkpoints
        else:
            self.checkpoints = None

        if self.forward and self.reverse:
            self.compounds = {**self.forward['compounds'], **self.reverse['compounds']}
            self.reactions = {**self.forward['reactions'], **self.reverse['reactions']}
            self.coreactants = {**self.forward['coreactants'], **self.reverse['coreactants']}
        elif self.forward:
            self.compounds = self.forward['compounds']
            self.reactions = self.forward['reactions']
            self.coreactants = self.forward['coreactants']
        elif self.reverse:
            self.compounds = self.reverse['compounds']
            self.reactions = self.reverse['reactions']
            self.coreactants = self.reverse['coreactants']

        if self.forward and not self.reverse and override_targets:
            self.forward['targets'] = override_targets
        elif self.reverse and not self.forward and override_starters:
            self.reverse['starters'] = override_starters
        elif override_starters or override_targets:
            raise ValueError("Cannot override starters or targets when both forward and reverse expansions are provided")

       
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
    def generation(self):
        if self.forward and self.reverse:
            return self.forward['generation'] + self.reverse['generation']
        elif self.forward:
            return self.forward['generation']
        elif self.reverse:
            return self.reverse['generation']

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
            'generation'
        ]
        half_expansion = {}

        with open(filepath, 'rb') as f:
            contents = pickle.load(f)

        for k in attributes:           
            if k == 'coreactants':
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
                reversed_reactions = dict([self._flip_reaction(rxn) for rxn in contents[k].values()])
                half_expansion[k] = reversed_reactions
            elif k == 'reactions':
                for key in contents[k].keys():
                    contents[k][key]['reversed'] = False # So that both fwd and rev half expansions have this field 
                half_expansion[k] = contents[k]
            else:
                half_expansion[k] = contents[k]

        if flip:
            half_expansion['targets'] = starters
            half_expansion['starters'] = targets
        else:
            half_expansion['starters'] = starters
            half_expansion['targets'] = targets

        return half_expansion

    def _flip_reaction(self, rxn: dict):
        flipped_rxn = {}
        flipped_rxn['_id'] = get_reaction_hash(rxn['Products'], rxn['Reactants'])
        flipped_rxn['Reactants'] = rxn['Products']
        flipped_rxn['Products'] = rxn['Reactants']
        flipped_rxn['Operators'] = rxn['Operators']
        flipped_rxn['SMILES_rxn'] = ' => '.join(rxn['SMILES_rxn'].split(' => ')[::-1])
        flipped_rxn['Operator_aligned_smarts'] = '>>'.join(rxn['Operator_aligned_smarts'].split('>>')[::-1])
        flipped_rxn['am_rxn'] = '>>'.join(rxn['am_rxn'].split('>>')[::-1])
        flipped_rxn['reversed'] = True
        return flipped_rxn['_id'], flipped_rxn
    
    def find_paths(self):
        '''
        Returns
        -------
        st_paths: list[list[list[str], list[str], list[str]]]
            Each element is a path. There are three subelements:
                - list of starter compound ids
                - list of target compound ids
                - list of reaction ids
            Everything is now in the "real-world" direction of synthesis
        '''
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

        st_paths = []
        for p in paths:
            s = set()
            t = set()
            rids = []
            for rid in p[1::2]:
                rids.append(rid)
                for _, cid in self.reactions[rid]['Reactants']:
                    if cid in self.starters:
                        s.add(cid)
                    
                for _, cid in self.reactions[rid]['Products']:
                    if cid in self.targets:
                        t.add(cid)

            st_paths.append([list(s), list(t), rids])
                    
        return st_paths    

    def _find_paths(self, half_expansion: dict, retro: bool):
        '''
        
        Args
        ----
        retro:bool
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
            paths += list(nx.all_simple_paths(DG, source=sr, target=sinks, cutoff=half_expansion['generation'] * 2)) # Search up to gens x 2 (bc bipartite)

        if flip:
            paths = [elt[::-1] for elt in paths]

        return paths
    
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
            
            # Neither starter, (mass source) nor X coreactant (non-mass-contributing source)
            non_sources = [
                c_id for _, c_id in rxn["Reactants"]
                if c_id[0] == 'C' and c_id not in self.starters
            ]

            mass_sources = [
                c_id for _, c_id in rxn["Reactants"]
                if c_id in self.starters
            ]

            if len(non_sources) == 1: # Other requirements for reaction are all sources
                edge_list.append((non_sources[0], i))
            elif len(non_sources) == 0: # Only sources required for reaction
                for ms in mass_sources:
                    edge_list.append((ms, i))
            else: # Currently don't support "extended branching" in synthesis paths
                pass

            for _, c_id in rxn["Products"]:

                # Don't outlink to non-mass-carrying coproducts
                # in order to approximately conserve mass along
                # synthesis paths
                if c_id[0] == 'X':
                     continue
                
                edge_list.append((i, c_id))
        
        # Create Graph
        DG = nx.DiGraph()  
        DG.add_nodes_from(cpd_node_list, bipartite=0)
        DG.add_nodes_from(rxn_node_list, bipartite=1)
        DG.add_edges_from(edge_list)
        
        return DG
    
    def prune(self, paths: list[list[list[str], list[str], list[str]]]):
        '''
        Prune compounds and reactions to save time
        e.g., with thermo stuff
        '''
        pruned_rxns = {}
        pruned_cpds = {}

        for _, _, rids in paths:
            for rid in rids:
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
    reactants: list[tuple], products: list[tuple]
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
    
    def __init__(self, study: pathlib.Path, known: pathlib.Path) -> None:
        self.study = study
        self.predicted_reactions = study / "predicted_reactions.parquet"
        self.paths = study / "found_paths.parquet"
        self.known_reactions = known / "known_reactions.parquet"
        self.enzymes = known / "known_enzymes.parquet"

        sts = pl.scan_parquet(self.paths).select(
            pl.col("starter_ids"),
            pl.col("target_ids"),
            pl.col("starters"),
            pl.col("targets")
        ).explode(
            [
                "starter_ids",
                "target_ids",
                "starters",
                "targets"
            ]
        ).unique()

        starters = sts.select(
            pl.col("starter_ids"),
            pl.col("starters")
        ).unique().collect()

        targets = sts.select(
            pl.col("target_ids"),
            pl.col("targets")
        ).unique().collect()

        self.starters = tuple(starters["starters"])
        self.targets = tuple(targets["targets"])
        self.starter_ids = tuple(starters["starter_ids"])
        self.target_ids = tuple(targets["target_ids"])

    def get_paths(
            self,
            starters: Iterable[str] = None, # TODO: switch to using ids
            targets: Iterable[str] = None,
            sort_by: str | None = None,
            descending: bool = True,
            lower_bounds: dict[str, float] | None = None,
            upper_bounds: dict[str, float] | None = None,
            filter_by_enzymes: dict[str, Iterable[str]] | None = None,
            top_k: int = None
        ) -> dict[str, pl.DataFrame]:
        '''
        Get valid paths with options to filter & sort

        Args
        -------
        starters: Iterable[str]
            Names of desired starter molecules
        targets: Iterable[str]
            Names of desired target molecules
        sort_by: str
            Field to sort by, can be one of:
                - "mdf"
                - "mean_max_rxn_sim"
                - "mean_mean_rxn_sim"
                - "min_max_rxn_sim"
                - "min_mean_rxn_sim"
                - "feasibility_frac"
        descending: bool = True
            If true sorts in descending order, otherwise ascending
        lower_bounds: dict[str, float]
            Maps field name -> lower bound value.
            Same fields as in sort_by.
        upper_bounds: dict[str, float]
            Maps field name -> upper bound value.
            Same fields as in sort_by.
        filter_by_enzymes: dict[str, Iterable]
            Maps field name -> list of desired values.
            Field may be one of:
                - "id"
                - "sequence"
                - "existence"
                - "reviewed"
                - "ec"
                - "organism"
                - "name"
        top_k: int
            If provided, returns top k paths sorted by sort_by field
        
        Returns
        ---------
        batch: dict[str, pl.DataFrame]
            {
                "paths": pl.DataFrame of paths,
                "predicted_reactions": pl.DataFrame of predicted reactions,
                "known_reactions": pl.DataFrame of known reactions,
                "enzymes": pl.DataFrame of enzymes
            }
        '''
        if top_k and not sort_by:
            raise ValueError("If top_k is provided, sort_by must also be provided")
        
        _path_schema = pl.read_parquet_schema(self.paths)
        _enzyme_schema = pl.read_parquet_schema(self.enzymes)

        if sort_by and sort_by not in _path_schema:
            raise ValueError(f"Invalid sort_by field: {sort_by}. Must be one of {_path_schema}")
        
        if lower_bounds and not all([f in _path_schema for f in lower_bounds]):
            raise ValueError(f"Invalid lower_bounds field(s): {lower_bounds}. Must be one of {_path_schema}")
        
        if upper_bounds and not all([f in _path_schema for f in upper_bounds]):
            raise ValueError(f"Invalid upper_bounds field(s): {upper_bounds}. Must be one of {_path_schema}")
        
        if filter_by_enzymes and not all([f in _enzyme_schema for f in filter_by_enzymes]):
            raise ValueError(f"Invalid filter_by_enzymes field(s): {filter_by_enzymes}. Must be one of {_enzyme_schema}")
        
        paths_lf = pl.scan_parquet(self.paths)

        if starters:
            paths_lf = paths_lf.filter(
                pl.col("starters").list.eval(pl.element().is_in(starters)).list.all()
            )
        
        if targets:
            paths_lf = paths_lf.filter(
                pl.col("targets").list.eval(pl.element().is_in(targets)).list.all()
            )

        if lower_bounds:
            for f, v in lower_bounds.items():
                paths_lf = paths_lf.filter(pl.col(f) >= v)
        
        if upper_bounds:
            for f, v in upper_bounds.items():
                paths_lf = paths_lf.filter(pl.col(f) <= v)

        if sort_by:
            paths_lf = paths_lf.sort(sort_by, descending=descending)

        paths_df = paths_lf.slice(0, top_k).collect()

        prids = set(paths_df['reactions'].explode()) # Get all unique reaction ids in paths
        prxns_df = pl.scan_parquet(self.predicted_reactions).filter(pl.col("id").is_in(prids)).collect()
        prxns_df = prxns_df.with_columns(
            pl.col("id").map_elements(lambda x: str(self.study / "svgs" / "predicted" / f"{x}.svg"), return_dtype=pl.String).alias("image"),
        )

        krids = set(prxns_df['analogue_ids'].explode()) # Get all unique known reaction ids in predicted reactions
        krxns_df = pl.scan_parquet(self.known_reactions).filter(pl.col("id").is_in(krids)).collect()
        krxns_df = krxns_df.with_columns(
            pl.col("id").map_elements(lambda x: str(self.study / "svgs" / "known" / f"{x}.svg"), return_dtype=pl.String).alias("image"),
        )

        enz_ids = set(krxns_df['enzymes'].explode()) # Get all unique enzyme ids in known reactions
        enz_lf = pl.scan_parquet(self.enzymes).filter(pl.col("id").is_in(enz_ids))

        if filter_by_enzymes:
            for f, v in filter_by_enzymes.items():
                enz_lf = enz_lf.filter(pl.col(f).is_in(v))

        enz_df = enz_lf.collect()

        return {"paths": paths_df, "predicted_reactions": prxns_df, "known_reactions": krxns_df, "enzymes": enz_df}

    # TODO: Reimplement w/ polars
    # def get_path_with_id(self, pid: str) -> Path:
    #     '''
    #     Get path with id pid
    #     '''
    #     # Required reaction ids for this path id
    #     req_prids = set()
    #     req_krids = set()
    #     for prid in self.paths[pid]['reactions']:
    #         req_prids.add(prid)
    #         for krid in self.predicted_reactions[prid]['analogues']:
    #             req_krids.add(krid)

    #     # Get known reactions
    #     krs = {}
    #     for krid in req_krids:
    #         kr = KnownReaction.from_dict(self.known_reactions[krid])
    #         krs[krid] = kr

    #     # Get predicted reactions
    #     prs = {}
    #     for prid in req_prids:
    #         pr = PredictedReaction.from_dict(self.predicted_reactions[prid], krs)
    #         prs[prid] = pr

    #     p = Path.from_dict(self.paths[pid], prs) # Construct path

    #     if not p.valid:
    #         print("Warning, this path may not have a complete set of known analogues & enzymes!")

    #     return p
    
    # def _validate_custom_agg(self, sort_by:dict):
    #     req_kwargs = {'kr_agg': str, 'pr_agg': str, 'k': int}
    #     for k, t in req_kwargs.items():
    #         if k not in sort_by:
    #             raise ValueError(f"Missing required argument {k}")
    #         if type(sort_by[k]) is not t:
    #             raise ValueError(f"Invalid type {sort_by[k]} for argument {k}. Must be {t}")

if __name__ == '__main__':
    from hydra import compose, initialize
    
    with initialize(version_base=None, config_path="conf/filepaths"):
        filepaths = compose(config_name="filepaths")
    
    forward = filepaths['raw_expansions'] / "2_steps_pivalic_acid_to_bottle_targets_24_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_True.pk"
    expansion = Expansion(forward=forward)
    paths = expansion.find_paths()

    forward = filepaths['raw_expansions'] / "2_steps_ccm_aa_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False_aplusb_True.pk"
    reverse = filepaths['raw_expansions'] / "2_steps_bottle_targets_24_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False_aplusb_False.pk"
    imt_reverses = load_json(filepaths['rules'] / "jnimt_reverses.json")
    expansion = Expansion(forward=forward, reverse=reverse, operator_reverses=imt_reverses)
    paths = expansion.find_paths()
    print()
