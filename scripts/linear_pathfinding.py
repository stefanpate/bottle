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
from networkx.exception import NodeNotFound

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

@hydra.main(version_base=None, config_path="../conf", config_name="linear_pathfinding")
def main(cfg: DictConfig):

    # Load data for reaction network
    cpds = pl.read_parquet(Path(cfg.expansion_extract) / cfg.cpds)
    helpers = cpds.filter(pl.col('type') == 'helper')['id'].to_list()
    sources = cpds.filter(pl.col('type') == 'source')['id'].to_list()
    targets = cpds.filter(pl.col('type') == 'target')['id'].to_list()
    cid2name = dict(zip(cpds['id'].to_list(), cpds['name'].to_list()))

    am_rxns = pl.read_parquet(Path(cfg.expansion_extract) / cfg.am_rxns)

    G = ReactionNetwork() # init
    logger.info("Adding reactions to network...")
    tic = perf_counter()
    for row in am_rxns.iter_rows(named=True):
        try:
            G.add_reaction(row['am_smarts'], rxn_type='predicted')
        except:
            logger.info(f"Failed to add reaction: {row['am_smarts']}")
            continue
    
    toc = perf_counter()
    logger.info(f"Added {len(am_rxns)} reactions in {toc - tic:.2f} seconds.")

    targets = [t for t in targets if G.has_node(t)]
    if len(targets) == 0:
        logger.info("No valid targets in network, exiting...")
        return
    
    logger.info("Setting helpers...")
    G.set_helpers(ids=helpers)

    
    logger.info("Finding paths...")
    paths = []
    tic = perf_counter()
    for sid in sources:
        if not G.has_node(sid):
            logger.info(f"Source {sid} not in network, skipping...")
            continue
        
        paths += nx.all_simple_edge_paths(
            G=G,
            source=sid,
            target=targets,
            cutoff=cfg.max_depth
        )

    toc = perf_counter()
    logger.info(f"Found {len(paths)} paths in {toc - tic:.2f} seconds.")
    
    logger.info("Saving results...")

    # Check for existing paths and reactions
    existing_reactions = check_existing("predicted_reactions.parquet", predicted_reactions_schema)
    existing_paths = check_existing("paths.parquet", paths_schema)
    existing_path_stats = check_existing("path_stats.parquet", path_stats_schema)

    # Lookups to generate new entries
    smarts_lookup = {}
    rxn_type_lookup = {}
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