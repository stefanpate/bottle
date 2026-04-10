import hydra
from omegaconf import DictConfig, OmegaConf
import polars as pl
from src.network import ReactionNetwork, de_am
from src.schemas import (
    predicted_reactions_schema,
    paths_schema,
    path_stats_schema,
    expansion_reactions_schema,
    compounds_schema,
    compound_type,
)
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
    starters, target_ids = [], []
    for generation, (rct, pdt, rid) in enumerate(path):
        rids.append(rid)
        gens.append(generation)
        main_pdt_ids.append(pdt)

        # Can have multiple starters/target_ids per path
        if rct in source_ids:
            starters.append(rct)
        if pdt in target_ids:
            target_ids.append(pdt)
            
    path_id = hash_path(list(zip(gens, rids)))

    return {
        "path_id": path_id,
        "rxn_ids": rids,
        "main_pdt_ids": main_pdt_ids,
        "generations": gens,
        "starters": starters,
        "target_ids": target_ids,
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
def construct_network(am_rxns: pl.DataFrame, source_ids: list[str], helper_ids: list[str], source_augmented_pnmc_lb: float, n_proc: int | None = None) -> ReactionNetwork:
    G = ReactionNetwork() # init
    logger.info("Adding reactions to network...")
    G.batch_add_reactions(
        am_rxns=am_rxns['am_smarts'].to_list(),
        rxn_types=['predicted'] * am_rxns.height,
        n_proc=n_proc,
    )
    logger.info(f"Added {len(am_rxns)} reactions")
    
    logger.info("Setting source_ids & helper_ids...")
    G.set_helpers(ids=helper_ids)
    G.set_sources(ids=source_ids)

    logger.info("Pruning network...")
    G.prune(augmented_pnmc_lb=source_augmented_pnmc_lb) # Connect compound to compound if it contributes at least 1.0 PNMC w/ source_ids+helper_ids
    return G

@timeit
def find_half_paths(G: nx.MultiDiGraph, source_ids: list[str], target_ids: list[str], max_depth: int):
    target_ids = [t for t in target_ids if G.has_node(t)]
    if len(target_ids) == 0:
        logger.info("No valid target_ids in network, exiting...")
        return []
    
    paths = []
    for sid in source_ids:
        if not G.has_node(sid):
            logger.info(f"Source {sid} not in network, skipping...")
            continue
        
        paths += nx.all_simple_edge_paths(
            G=G,
            source=sid,
            target=target_ids,
            cutoff=max_depth
        )

    return paths

@timeit
def find_combo_paths(graphs: list[nx.MultiDiGraph], source_ids: list[str], checkpoint_ids: list[str], target_ids: list[str], max_depths: dict[str, int]):
    if len(graphs) != len(max_depths):
        raise ValueError("Length of graphs and max_depths must be the same.")
    
    # Remove source_ids from sinks otherwise find paths is ill defined and nx fails
    # Union w/ target_ids to catch target reached in half expansion
    _checkpoint_ids = [c for c in set(checkpoint_ids) | set(target_ids) if c not in source_ids]
    fwd_paths = find_half_paths(graphs[0], source_ids, _checkpoint_ids, max_depths['forward'])

    # Remove target_ids from checkpoint_ids otherwise find paths is ill defined and nx fails
    # Union w/ source_ids to catch source reached in half expansion
    _checkpoint_ids = [c for c in set(checkpoint_ids) | set(source_ids) if c not in target_ids]
    rev_paths = find_half_paths(graphs[-1], _checkpoint_ids, target_ids, max_depths['retro'])
    
    all_paths = []
    fwd_by_checkpoint = defaultdict(list)
    for fp in fwd_paths:
        if fp[-1][1] in target_ids: # Final product is a target, done
            all_paths.append(fp)
        else:
            fwd_by_checkpoint[fp[-1][1]].append(fp) # Save to stitch with reverse paths

    rev_by_checkpoint = defaultdict(list)
    for rp in rev_paths:
        if rp[0][0] in source_ids: # Initial reactant is a source, done
            all_paths.append(rp)
        else:
            rev_by_checkpoint[rp[0][0]].append(rp) # Save to stitch with forward paths

    # Stitch together at checkpoint_ids
    for ckpt in fwd_by_checkpoint.keys() & rev_by_checkpoint.keys():
        for fp, rp in product(fwd_by_checkpoint[ckpt], rev_by_checkpoint[ckpt]):
            all_paths.append(fp + rp)

    return all_paths            


@hydra.main(version_base=None, config_path="../conf", config_name="linear_pathfinding")
def main(cfg: DictConfig):

    # Load data for reaction network
    if cfg.forward_expansion:
        fwd_rxns = pl.read_parquet(
            Path(cfg.fwd_dir) / cfg.am_rxns, schema=expansion_reactions_schema
        ).with_columns(
            pl.lit('forward').alias('half_expansion')
        )
        fwd_cpds = pl.read_parquet(Path(cfg.fwd_dir) / cfg.cpds, schema=compounds_schema)
    else:
        fwd_rxns = pl.DataFrame(
            schema=expansion_reactions_schema
        ).with_columns(
            pl.lit('forward').alias('half_expansion')
        )
        fwd_cpds = pl.DataFrame(schema=compounds_schema)

    if cfg.retro_expansion:
        retro_rxns = pl.read_parquet(
            Path(cfg.retro_dir) / cfg.am_rxns, schema=expansion_reactions_schema
        ).with_columns(
            pl.lit('retro').alias('half_expansion')
        )
        retro_cpds = pl.read_parquet(Path(cfg.retro_dir) / cfg.cpds, schema=compounds_schema)
    else:
        retro_rxns = pl.DataFrame(
            schema=expansion_reactions_schema
        ).with_columns(
            pl.lit('retro').alias('half_expansion')
        )
        retro_cpds = pl.DataFrame(schema=compounds_schema)

    if cfg.forward_expansion and cfg.retro_expansion:
        checkpoint_ids = fwd_cpds.join(
            other=retro_cpds,
            on='id',
            how='inner'
        )['id'].to_list()
    else:
        checkpoint_ids = []

    am_rxns = pl.concat([fwd_rxns, retro_rxns]).unique()
    cpds = pl.concat([fwd_cpds, retro_cpds]).with_columns(
        pl.when(
            (pl.col('type') == 'intermediate') & (pl.col('id').is_in(checkpoint_ids))
        ).then(
            pl.lit('checkpoint')
        ).otherwise(
            pl.col('type')
        ).alias('type').cast(compound_type)
    )

    sources = cpds.filter(pl.col('type') == 'source')
    targets = cpds.filter(pl.col('type') == 'target')
    helper_ids = cpds.filter(pl.col('type') == 'helper')['id'].to_list()
    source_ids = sources['id'].to_list()
    target_ids = targets['id'].to_list()
    # checkpoint_ids = cpds.filter(pl.col('type') == 'checkpoint')['id'].to_list()
    sid2name = dict(zip(sources['id'], sources['name']))
    tid2name = dict(zip(targets['id'], targets['name']))

    if cfg.forward_expansion and cfg.retro_expansion:
        mode = 'combo'
    elif cfg.forward_expansion:
        mode = 'forward'
    elif cfg.retro_expansion:
        mode = 'retro'

    generations = {}
    if cfg.forward_expansion:
        with open(Path(cfg.fwd_dir) / cfg.expansion_config, 'r') as f:
            fwd_cfg = OmegaConf.load(f)
            generations['forward'] = fwd_cfg.generations

    if cfg.retro_expansion:
        with open(Path(cfg.retro_dir) / cfg.expansion_config, 'r') as f:
            retro_cfg = OmegaConf.load(f)
            generations['retro'] = retro_cfg.generations

    logger.info(f"Mode: {mode}")
    if mode == 'forward' or mode == 'retro':
        G = construct_network(am_rxns, source_ids, helper_ids, cfg.source_augmented_pnmc_lb, cfg.processes)
        logger.info("Finding paths...")
        paths = find_half_paths(G, source_ids, target_ids, generations[mode])
        Gs = [G]
    elif mode == 'combo':
        F = construct_network(am_rxns.filter(pl.col('half_expansion') == 'forward'), source_ids, helper_ids, cfg.source_augmented_pnmc_lb, cfg.processes) # Setting addtl mass source_ids but
        R = construct_network(am_rxns.filter(pl.col('half_expansion') == 'retro'), source_ids, helper_ids, cfg.source_augmented_pnmc_lb, cfg.processes) # Setting addtl mass source_ids but still pathfinding from checkpoint_ids to target_ids
        Gs = [F, R]
        logger.info("Finding paths...")
        paths = find_combo_paths(Gs, source_ids, checkpoint_ids, target_ids, generations)

    logger.info(f"Found {len(paths)} paths")
    
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

    # Collect paths and path stats
    new_paths, new_path_stats = [], []
    new_rxn_ids = set()
    n_existing = 0
    for path in paths:
        path_entry = path_to_path_entry(path, source_ids, target_ids)
        
        if path_entry['path_id'] in existing_paths['path_id'].to_list():
            n_existing += 1
            continue

        # One entry per path
        new_path_stats.append(
            [
                path_entry['path_id'],
                [sid2name.get(sid) for sid in path_entry['starters']],
                [tid2name.get(tid) for tid in path_entry['target_ids']],
                None, # dg_opt
                None, # dg_err
                path_entry['starters'], # starter_ids
                path_entry['target_ids'], # target_ids
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

    logger.info(f"{n_existing} / {len(paths)} paths were already stored")
    logger.info(f"Total # new paths found: {len(new_path_stats)}")
    
    # Merge identical new paths
    to_merge = [
        'starters',
        'target_ids',
        'starter_ids',
        'target_ids',
    ]
    aggs = []
    for col in path_stats_schema.names():
        if col == "id":
            continue
        
        if col in to_merge:
            aggs.append(pl.col(col).flatten().unique().alias(col))
        else:
            aggs.append(pl.col(col).first().alias(col))

    new_path_stats = pl.DataFrame(
        new_path_stats,
        schema=path_stats_schema,
        orient='row'
    ).group_by("id").agg(aggs)

    new_paths = pl.DataFrame(
        new_paths,
        schema=paths_schema,
        orient='row'
    ).unique()

    # Concat paths, path_stats
    
    paths = pl.concat([existing_paths, new_paths])
    path_stats = pl.concat([existing_path_stats, new_path_stats])

    # Collect new unique predicted reactions
    new_reactions = []
    for rxn_id in new_rxn_ids:

        if rxn_type_lookup[rxn_id] == 'predicted' and rxn_id not in existing_reactions['id']:
            am_rxn = smarts_lookup[rxn_id]
            tmp = am_rxns.filter(pl.col('am_smarts') == am_rxn)
            rules = tmp['rule_name'].to_list()
            rule_sets = tmp['rule_set'].to_list()
            templates = tmp['rule_template'].to_list()
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
                    rules, # rule names
                    templates, # templates
                    rule_sets, # rule_sets
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
    logger.info(f"Saving results to {cfg.casp_study}")
    paths.write_parquet("paths.parquet")
    path_stats.write_parquet("path_stats.parquet")
    reactions.write_parquet("predicted_reactions.parquet")
    cpds.write_parquet("compounds.parquet")

if __name__ == "__main__":
    main()