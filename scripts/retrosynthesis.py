import hydra
from omegaconf import DictConfig
import polars as pl
from src.network import ReactionNetwork
from src.schemas import predicted_reactions_schema, paths_schema, path_stats_schema
from pathlib import Path
from time import perf_counter
import logging

logger  = logging.getLogger(__name__)

def check_existing(fn: str, schema: pl.Schema):
    if (dir / fn).exists():
        existing = pl.scan_parquet(
            dir / fn
        ).select(
            pl.col("id")
        ).collect()
    else:
        existing = pl.DataFrame(schema=schema)
        existing.write_parquet(dir / fn)
    return existing

@hydra.main(version_base=None, config_path="../configs", config_name="retrosynthesis")
def main(cfg: DictConfig):

    # Load data for reaction network

    G = ReactionNetwork.from_json(Path(cfg.known_reaction_network))

    default_sources = pl.read_csv(Path(cfg.default_sources))['smiles'].to_list()

    with open(Path(cfg.expansion_extract) / cfg.sources, 'r') as f:
        sources = [line.strip() for line in f.readlines()]

    with open(Path(cfg.expansion_extract) / cfg.targets, 'r') as f:
        targets = [line.strip() for line in f.readlines()]

    am_rxns = pl.read_parquet(Path(cfg.expansion_extract) / cfg.am_rxns)

    # TODO: Add operators as field
    logger.info("Adding reactions to network...")
    tic = perf_counter()
    for row in am_rxns.iter_rows(named=True):
        try:
            G.add_reaction(row['am_smarts'])
        except:
            logger.info(f"Failed to add reaction: {row['am_smarts']}")
            continue
    toc = perf_counter()
    logger.info(f"Added {len(am_rxns)} reactions in {toc - tic:.2f} seconds.")

    rxn_lookup = {k: d['am_smarts'] for _, _, k, d in G.edges(keys=True, data=True)}

    logger.info("Setting sources...")
    G.set_sources(smiles=default_sources)
    G.set_sources(smiles=sources)

    target_ids = [G.get_nodes_by_prop('smiles', t)[0] for t in targets]

    logger.info("Enumerating synthetic trees...")
    trees = []
    for tid in target_ids:
        tic = perf_counter()
        trees.append(
            G.enumerate_synthetic_trees(
                target=tid,
                max_depth=cfg.max_depth,
                max_leaves=cfg.max_leaves,
                tot_rnmc_lb=cfg.tot_rnmc_lb
            )
        )

        toc = perf_counter()
        logger.info(f"Tree enumeration for target {tid} took {toc - tic:.2f} seconds.")
    
    logger.info("Saving results...")

    # Check for existing paths and reactions
    existing_reactions = check_existing("predicted_reactions.parquet", predicted_reactions_schema)
    existing_paths = check_existing("paths.parquet", paths_schema)
    existing_path_stats = check_existing("path_stats.parquet", path_stats_schema)
    

    with open(f"{cfg.expansion}_synthetic_trees.json", 'w') as f:
        json.dump(target_to_trees, f)

if __name__ == "__main__":
    main()