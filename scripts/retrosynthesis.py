import hydra
from omegaconf import DictConfig, OmegaConf
import polars as pl
from src.network import ReactionNetwork, SyntheticTree, de_am
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

logger  = logging.getLogger(__name__)

def check_existing(fn: str, schema: pl.Schema):
    if (Path.cwd() / fn).exists():
        existing = pl.read_parquet(Path.cwd() / fn)
    else:
        existing = pl.DataFrame(schema=schema)
    return existing

def tree_to_path_entry(tree: SyntheticTree):
    rids, gens, main_pdt_ids = [], [], []
    for i, gen in enumerate(tree.generations):
        for main_pdt_id, rxn_id in gen.items():
            
            if rxn_id is None:
                continue
            
            rids.append(rxn_id)
            gens.append(i)
            main_pdt_ids.append(main_pdt_id)
    path_id = hash_path(list(zip(gens, rids)))
    target = tree.root
    starters = [cid for cid, _ in tree.leaves]

    return {
        "path_id": path_id,
        "rxn_ids": rids,
        "main_pdt_ids": main_pdt_ids,
        "generations": gens,
        "starters": starters,
        "targets": [target],
    }

@hydra.main(version_base=None, config_path="../conf", config_name="retrosynthesis")
def main(cfg: DictConfig):

    # Load data for reaction network
    G = ReactionNetwork.from_json(Path(cfg.known_reaction_network))
    default_sources = pl.read_csv(Path(cfg.default_sources))['smiles'].to_list()

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

    am_rxns = pl.concat([fwd_rxns, retro_rxns]).unique()
    cpds = pl.concat([fwd_cpds, retro_cpds]).unique()

    helpers = cpds.filter(pl.col('type') == 'helper')['id'].to_list()
    sources = cpds.filter(pl.col('type') == 'source')['id'].to_list()
    targets = cpds.filter(pl.col('type') == 'target')['id'].to_list()
    cid2name = dict(zip(cpds['id'].to_list(), cpds['name'].to_list()))

    generations = {}
    if cfg.forward_expansion:
        with open(Path(cfg.fwd_dir) / cfg.expansion_config, 'r') as f:
            fwd_cfg = OmegaConf.load(f)
            generations['forward'] = fwd_cfg.generations

    if cfg.retro_expansion:
        with open(Path(cfg.retro_dir) / cfg.expansion_config, 'r') as f:
            retro_cfg = OmegaConf.load(f)
            generations['retro'] = retro_cfg.generations


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

    logger.info("Setting sources...")
    G.set_sources(ids=sources)
    G.set_helpers(smiles=default_sources)
    G.set_helpers(ids=helpers)

    logger.info("Enumerating synthetic trees...")
    trees = []
    for tid in targets:
        if not G.has_node(tid):
            logger.info(f"Target {tid} not in reaction network, skipping...")
            continue

        tic = perf_counter()
        trees += G.enumerate_synthetic_trees(
                target=tid,
                max_depth=sum(generations.values()),
                max_leaves=cfg.max_leaves,
                tot_rnmc_lb=cfg.tot_rnmc_lb
            )

        toc = perf_counter()
        logger.info(f"Tree enumeration for target {tid} took {toc - tic:.2f} seconds.")
    

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

    # Collect paths and path stats

    new_paths, new_path_stats = [], []
    new_rxn_ids = set()
    n_existing = 0
    for tree in trees:
        path_entry = tree_to_path_entry(tree)
        if path_entry['path_id'] in existing_paths['path_id'].to_list():
            n_existing += 1
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

    logger.info(f"{n_existing} / {len(trees)} paths were already stored")
    logger.info(f"Total # new paths found: {len(new_path_stats)}")

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

if __name__ == "__main__":
    main()