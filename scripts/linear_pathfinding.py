import hydra
from omegaconf import DictConfig
import polars as pl
from src.network import ReactionNetwork, de_am
from src.schemas import predicted_reactions_schema, paths_schema, path_stats_schema
from src.post_processing import hash_path
from pathlib import Path
from time import perf_counter
import logging
import networkx as nx
from collections import defaultdict
from itertools import product
from functools import wraps

logger  = logging.getLogger(__name__)

def check_existing(fn: str, schema: pl.Schema):
    if (Path.cwd() / fn).exists():
        existing = pl.read_parquet(Path.cwd() / fn)
    else:
        existing = pl.DataFrame(schema=schema)
    return existing

def path_to_path_entry(path: list[tuple[str, str, str]], source_ids: list[str], target_ids: list[str]):
    rids, gens, main_pdt_ids = [], [], []
    starters, targets = [], []
    for generation, (rct, pdt, rid) in enumerate(path):
        rids.append(rid)
        gens.append(generation)
        main_pdt_ids.append(pdt)

        # Can have multiple starters/targets per path
        if rct in source_ids:
            starters.append(rct)
        if pdt in target_ids:
            targets.append(pdt)
            
    path_id = hash_path(list(zip(gens, rids)))

    return {
        "path_id": path_id,
        "rxn_ids": rids,
        "main_pdt_ids": main_pdt_ids,
        "generations": gens,
        "starters": starters,
        "targets": targets,
    }


def timeit(fcn):
    @wraps(fcn)
    def wrapper(*args, **kwargs):
        tic = perf_counter()
        result = fcn(*args, **kwargs)
        toc = perf_counter()
        logger.info(f"Function {fcn.__name__} took {toc - tic:.2f} seconds.")
        return result
    return wrapper

@timeit
def construct_network(am_rxns: pl.DataFrame, helpers: list[str]) -> ReactionNetwork:
    G = ReactionNetwork() # init
    logger.info("Adding reactions to network...")
    for row in am_rxns.iter_rows(named=True):
        try:
            G.add_reaction(row['am_smarts'], rxn_type='predicted')
        except:
            logger.info(f"Failed to add reaction: {row['am_smarts']}")
            continue
    
    logger.info(f"Added {len(am_rxns)} reactions")
    
    logger.info("Setting helpers...")
    G.set_helpers(ids=helpers)

    logger.info("Pruning network...")
    G.prune()
    return G

@timeit
def find_half_paths(G: nx.MultiDiGraph, sources: list[str], targets: list[str], max_depth: int):
    targets = [t for t in targets if G.has_node(t)]
    if len(targets) == 0:
        logger.info("No valid targets in network, exiting...")
        return []
    
    paths = []
    for sid in sources:
        if not G.has_node(sid):
            logger.info(f"Source {sid} not in network, skipping...")
            continue
        
        paths += nx.all_simple_edge_paths(
            G=G,
            source=sid,
            target=targets,
            cutoff=max_depth
        )

    return paths

@timeit
def find_combo_paths(graphs: list[nx.MultiDiGraph], sources: list[str], checkpoints: list[str], targets: list[str], max_depths: dict[str, int]):
    if len(graphs) != len(max_depths):
        raise ValueError("Length of graphs and max_depths must be the same.")
    
    # Remove sources from sinks otherwise find paths is ill defined and nx fails
    # Union w/ targets to catch target reached in half expansion
    _checkpoints = [c for c in set(checkpoints) | set(targets) if c not in sources]
    fwd_paths = find_half_paths(graphs[0], sources, _checkpoints, max_depths['forward'])

    # Remove targets from checkpoints otherwise find paths is ill defined and nx fails
    # Union w/ sources to catch source reached in half expansion
    _checkpoints = [c for c in set(checkpoints) | set(sources) if c not in targets]
    rev_paths = find_half_paths(graphs[-1], _checkpoints, targets, max_depths['retro'])
    
    all_paths = []
    fwd_by_checkpoint = defaultdict(list)
    for fp in fwd_paths:
        if fp[-1][-1] in targets: # Final product is a target, done
            all_paths.append(fp)
        else:
            fwd_by_checkpoint[fp[-1][-1]].append(fp) # Save to stitch with reverse paths

    rev_by_checkpoint = defaultdict(list)
    for rp in rev_paths:
        if rp[0][0] in sources: # Initial reactant is a source, done
            all_paths.append(rp)
        else:
            rev_by_checkpoint[rp[0][0]].append(rp) # Save to stitch with forward paths

    # Stitch together at checkpoints
    for ckpt in fwd_by_checkpoint.keys() & rev_by_checkpoint.keys():
        for fp, rp in product(fwd_by_checkpoint[ckpt], rev_by_checkpoint[ckpt]):
            all_paths.append(fp + rp)

    return all_paths            


@hydra.main(version_base=None, config_path="../conf", config_name="linear_pathfinding")
def main(cfg: DictConfig):

    # Load data for reaction network
    cpds = pl.read_parquet(Path(cfg.expansion_extract) / cfg.cpds)
    helpers = cpds.filter(pl.col('type') == 'helper')['id'].to_list()
    sources = cpds.filter(pl.col('type') == 'source')['id'].to_list()
    targets = cpds.filter(pl.col('type') == 'target')['id'].to_list()
    checkpoints = cpds.filter(pl.col('type') == 'checkpoint')['id'].to_list()
    cid2name = dict(zip(cpds['id'].to_list(), cpds['name'].to_list()))

    am_rxns = pl.read_parquet(Path(cfg.expansion_extract) / cfg.am_rxns)
    half_expansions = am_rxns['half_expansion'].unique().to_list()
    mode = half_expansions[0] if len(half_expansions) == 1 else 'combo'

    generations = pl.read_csv(Path(cfg.expansion_extract) / cfg.generations)
    generations = dict(zip(generations['half_expansion'].to_list(), generations['generation'].to_list()))

    logger.info(f"Mode: {mode}")
    if mode == 'forward' or mode == 'retro':
        G = construct_network(am_rxns, helpers)
        logger.info("Finding paths...")
        paths = find_half_paths(G, sources, targets, generations[mode])
        Gs = [G]
    elif mode == 'combo':
        F = construct_network(am_rxns.filter(pl.col('half_expansion') == 'forward'), helpers)
        R = construct_network(am_rxns.filter(pl.col('half_expansion') == 'retro'), helpers)
        Gs = [F, R]
        logger.info("Finding paths...")
        paths = find_combo_paths(Gs, sources, checkpoints, targets, generations)

    logger.info(f"Found {len(paths)} paths")
    
    logger.info("Saving results...")

    # Check for existing paths and reactions
    existing_reactions = check_existing("predicted_reactions.parquet", predicted_reactions_schema)
    existing_paths = check_existing("paths.parquet", paths_schema)
    existing_path_stats = check_existing("path_stats.parquet", path_stats_schema)

    # Lookups to generate new entries
    smarts_lookup = {}
    rxn_type_lookup = {}
    for G in Gs:
        for _, _, k, d in G.edges(keys=True, data=True):
            smarts_lookup[k] = d['am_smarts']
            rxn_type_lookup[k] = d['rxn_type']

    rules_lookup = dict(zip(am_rxns['am_smarts'], am_rxns['rules'].to_list()))

    # Collect paths and path stats
    new_paths, new_path_stats = [], []
    new_rxn_ids = set()
    for path in paths:
        path_entry = path_to_path_entry(path, sources, targets)
        
        if path_entry['path_id'] in existing_paths['path_id'].to_list():
            continue

        # One entry per path
        new_path_stats.append(
            [
                path_entry['path_id'],
                [cid2name.get(sid) for sid in path_entry['starters']],
                [cid2name.get(tid) for tid in path_entry['targets']],
                None, # dg_opt
                None, # dg_err
                path_entry['starters'], # starter_ids
                path_entry['targets'], # target_ids
                None, # mdf
                None, # mean_max_rxn_sim
                None, # mean_mean_rxn_sim
                None, # min_max_rxn_sim
                None, # min_mean_rxn_sim
                None, # feasibility_frac
            ]
        )

        # One entry per reaction
        for rxn_id, main_pdt_id, generation in zip(path_entry['rxn_ids'], path_entry['main_pdt_ids'], path_entry['generations']):
            new_paths.append(
                [
                    path_entry['path_id'],
                    rxn_id,
                    main_pdt_id,
                    rxn_type_lookup[rxn_id],
                    generation
                ]
            )
            new_rxn_ids.add(rxn_id)

    # Concat paths, path_stats
    
    paths = pl.concat([
        existing_paths,
        pl.DataFrame(
            new_paths,
            schema=paths_schema,
            orient='row'
        )
    ])

    path_stats = pl.concat([
        existing_path_stats,
        pl.DataFrame(
            new_path_stats,
            schema=path_stats_schema,
            orient='row'
        )
    ])

    # Collect new unique predicted reactions
    new_reactions = []
    for rxn_id in new_rxn_ids:

        if rxn_type_lookup[rxn_id] == 'predicted' and rxn_id not in existing_reactions['id']:
            am_rxn = smarts_lookup[rxn_id]
            rcts, pdts = de_am(am_rxn)
            de_am_rxn = f"{'.'.join(rcts)}>>{'.'.join(pdts)}"
            new_reactions.append(
                [
                    rxn_id,
                    de_am_rxn,
                    am_rxn,
                    None, # dxgb_label
                    None, # rxn_sims
                    None, # analogue_ids
                    rules_lookup[smarts_lookup[rxn_id]], # rules
                ]
            )

    
    # Concat and save predicted reactions
    reactions = pl.concat([
        existing_reactions,
        pl.DataFrame(
            new_reactions,
            schema=predicted_reactions_schema,
            orient='row'
        )
    ])

    # Save
    paths.write_parquet("paths.parquet")
    path_stats.write_parquet("path_stats.parquet")
    reactions.write_parquet("predicted_reactions.parquet")
    cpds.write_parquet("compounds.parquet")

if __name__ == "__main__":
    main()