import hydra
from omegaconf import DictConfig
from tqdm import tqdm
import polars as pl
from pathlib import Path
from src.schemas import path_stats_schema
from ergochemics.standardize import hash_compound
import numpy as np
from time import perf_counter
from collections import defaultdict
from logging import getLogger
from equilibrator_assets.local_compound_cache import LocalCompoundCache


def update_table(existing: pl.DataFrame, analyzed: pl.DataFrame, on: str = "id") -> pl.DataFrame:
    updated = existing.join(
        analyzed, on=on, how='left', suffix='_analyzed'
    )
    for col in existing.columns:
        analyzed_col = f"{col}_analyzed"
        if analyzed_col in updated.columns:
            updated = updated.with_columns(
                pl.col(col).fill_null(pl.col(analyzed_col)).alias(col)
            )
    updated = updated.select(existing.columns)
    return updated

# def proc_initializer(cfg: DictConfig, _mapped_rxns: pl.DataFrame) -> None:
#     print("", flush=True)
#     global 
    
logger = getLogger(__name__)

@hydra.main(version_base=None, config_path="../conf", config_name="analyze_thermo")
def main(cfg: DictConfig) -> None:

    if not Path("paths.parquet").exists() or not Path("predicted_reactions.parquet").exists() or not Path("path_stats.parquet").exists():
        return
    
    # Load semi-processed data

    cpds = pl.read_parquet("compounds.parquet")
    cid2name = dict(zip(cpds['id'].to_list(), cpds['name'].to_list()))

    path_stats_to_do = pl.scan_parquet("path_stats.parquet").filter(
        pl.col("mdf").is_null() | pl.col("dg_opt").is_null() | pl.col("dg_err").is_null()
    ).collect()

    paths_to_do = pl.scan_parquet("paths.parquet").filter(
        pl.col("path_id").is_in(path_stats_to_do["id"].to_list())
    ).select(
        pl.col("path_id"),
        pl.col("rxn_id"),
        pl.col("generation"),
    ).collect()
    
    pred_rxns_to_do = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("id").is_in(paths_to_do["rxn_id"].unique().to_list())
    ).select(
        pl.col("id"),
        pl.col("smarts"),
    ).collect()

    known_rxns_to_do = pl.scan_parquet(
        Path(cfg.filepaths.known_reactions)
    ).filter(
        pl.col("id").is_in(paths_to_do["rxn_id"].unique().to_list())
    ).select(
        pl.col("id"),
        pl.col("smarts"),
    ).collect()

    rxns_to_do = pl.concat((pred_rxns_to_do, known_rxns_to_do))

    # Get unique compounds
    unique_cpds = set()
    for smarts in rxns_to_do['smarts']:
        lhs, rhs = smarts.split(">>")
        for smi in lhs.split(".") + rhs.split("."):
            unique_cpds.add(smi)

    unique_cpds = pl.DataFrame(
        {
            "struct": list(unique_cpds),
            "coco_id": [hash_compound(smi) for smi in unique_cpds],
            "name": [cid2name.get(smi, "unknown") for smi in unique_cpds],
        }
    ).to_pandas()

    logger.info(f"{len(unique_cpds)} unique compounds / {len(rxns_to_do)} reactions / {len(path_stats_to_do)} paths to analyze")

    logger.info("Adding compounds to local cache...")
    lc = LocalCompoundCache()
    lc.load_cache(cfg.eq_cache)
    start = perf_counter()
    lc.add_compounds(
        unique_cpds,
        mol_format="smiles",
        bypass_chemaxon=True,
        save_empty_compounds=True,
    )
    end = perf_counter()
    logger.info(f"Added {len(unique_cpds)} compounds to local cache in {end - start:.2f} seconds")

    # # Add reaction-derived summary stats to paths
    # analyzed_path_stats = []
    # for row in tqdm(paths_to_do.iter_rows(named=True), total=len(paths_to_do), desc="Updating path stats"):
    #     # If rxn_sims (& analogue_ids) for a predicted rxn is empty, there are 
    #     # no analogues. When calculating summary stats, treat this as rxn_sims = [0.0]
    #     # instead of skipping it for example.

    #     rxns = req_rxns.filter(
    #         pl.col("id").is_in(row['rxn_id'])
    #     ).with_columns(
    #         pl.col("rxn_sims").list.mean().fill_null(0).alias("mean_rxn_sims"),
    #         pl.col("rxn_sims").list.max().fill_null(0).alias("max_rxn_sims"),
    #     )

    #     analyzed_path_stats.append(
    #         {
    #             "id": row['path_id'],
    #             "feasibility_frac": rxns["dxgb_label"].mean(),
    #             "mean_max_rxn_sim": rxns["max_rxn_sims"].mean(),
    #             "mean_mean_rxn_sim": rxns["mean_rxn_sims"].mean(),
    #             "min_max_rxn_sim": rxns["max_rxn_sims"].min(),
    #             "min_mean_rxn_sim": rxns["mean_rxn_sims"].min(),
    #         }
    #     )

    # analyzed_path_stats = pl.from_dicts(analyzed_path_stats, schema=path_stats_schema)
    # existing_path_stats = pl.read_parquet("path_stats.parquet")
    # updated_path_stats = update_table(existing_path_stats, analyzed_path_stats)
    # updated_path_stats.write_parquet("path_stats.parquet")

if __name__ == "__main__":
    main()