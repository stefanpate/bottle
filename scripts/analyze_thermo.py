import hydra
from omegaconf import DictConfig
from tqdm import tqdm
import polars as pl
from concurrent.futures import ProcessPoolExecutor
from hydra.utils import instantiate
from functools import partial
from itertools import chain
from ergochemics.mapping import get_reaction_center
from pathlib import Path
from rdkit import Chem
from src.schemas import path_stats_schema, predicted_reactions_schema
from src.chem_draw import draw_reaction
import numpy as np
from time import perf_counter
from typing import Any
from multiprocessing import get_context
from collections import defaultdict

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

def proc_initializer(cfg: DictConfig, _mapped_rxns: pl.DataFrame) -> None:
    def sma_rc_fp(sma, rc):
        mol = Chem.MolFromSmiles(sma.split('>>')[0])
        if mol is None:
            return None
        return _fingerprint(mol, rc)
    
    print("Initializing process for reaction mapping", flush=True)
    global dxgb, _fingerprint, mapped_rxns, rxn_reverses
    dxgb = instantiate(cfg.dxgb)
    mfper = instantiate(cfg.mfper)
    _fingerprint = partial(mfper.fingerprint, rc_dist_ub=cfg.rc_dist_ub)

    # Update mapped reactions with reaction center and fingerprints
    mapped_rxns = _mapped_rxns.with_columns(
        reaction_center=pl.col("am_smarts").map_elements(
            get_lhs_block_rc, return_dtype=pl.List(pl.Int64)
        )
    ).with_columns(
        pl.struct(["smarts", "reaction_center"]).map_elements(
            lambda x: sma_rc_fp(x["smarts"], x["reaction_center"]),
            return_dtype=pl.Object
        ).alias("mfp")
    )

def get_lhs_block_rc(am_smarts: str) -> list[int]:
    return get_reaction_center(am_smarts, mode="combined")[0]

def process_reaction(reaction: dict[str, Any]) -> dict[str, Any]:
    query_smarts = reaction["smarts"]
    query_am_smarts = reaction["am_smarts"]
    query_lhs_mol = Chem.MolFromSmiles(query_smarts.split('>>')[0])
    query_lhs_block_rc = get_lhs_block_rc(query_am_smarts)

    if query_lhs_mol is None:
        return reaction
    
    analogues = mapped_rxns.filter(pl.col("rule_id").is_in(reaction["rules"]))

    if analogues.is_empty():
        srt_sims = []
        srt_krids = []
    else:
        mfps = np.vstack(analogues["mfp"])
        query_mfp = _fingerprint(
            mol=query_lhs_mol,
            reaction_center=query_lhs_block_rc
        ).reshape(-1, 1)
        sims = np.matmul(mfps, query_mfp) / (np.linalg.norm(mfps, axis=1).reshape(-1, 1) * np.linalg.norm(query_mfp))
        sims = sims.reshape(-1,)
        srt_sims = sorted(sims, reverse=True)
        srt_idxs = np.argsort(sims)[::-1]
        srt_krids = analogues['rxn_id'][srt_idxs].to_list()

    is_feasible = dxgb.predict_label(query_smarts)
    reaction["dxgb_label"] = is_feasible
    reaction["rxn_sims"] = srt_sims
    reaction["analogue_ids"] = srt_krids
    return reaction

def load_mapped_reactions(rule_names: pl.Series, mappings_dir: str) -> pl.DataFrame:
    unique_rules = set(rule_names.explode())
    rules_by_set = defaultdict(set)
    for rule in unique_rules:
        ruleset_name, rule_id = rule.split(":")
        rules_by_set[ruleset_name].add(int(rule_id))
    
    mapped_rxns = []
    for ruleset_name, rule_ids in rules_by_set.items():
        df = pl.scan_parquet(
            Path(mappings_dir) / f"mapped_known_reactions_x_{ruleset_name}_rules.parquet"
        ).filter(
            pl.col("rule_id").is_in(rule_ids)
        ).with_columns(
            pl.col("rule_id").map_elements(lambda x: f"{ruleset_name}:{x}", return_dtype=pl.String).alias("rule_id")
        ).collect()
        mapped_rxns.append(df)

    mapped_rxns = pl.concat(mapped_rxns)
    return mapped_rxns


@hydra.main(version_base=None, config_path="../conf", config_name="analyze_structures")
def main(cfg: DictConfig) -> None:

    if not Path("paths.parquet").exists():
        return
    
    pred_rxns_to_do = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("dxgb_label").is_null() | pl.col("rxn_sims").is_null() | pl.col("analogue_ids").is_null()
    ).collect()

    if not len(pred_rxns_to_do) == 0:
        _mapped_rxns = load_mapped_reactions(rule_names=pred_rxns_to_do["rules"], mappings_dir=cfg.filepaths.rxn_x_rule_mapping)
        
        # Process reactions
        chunksize = max(1, int(len(pred_rxns_to_do) / cfg.processes))
        context = get_context("spawn") # Polars hates fork
        with ProcessPoolExecutor(max_workers=cfg.processes, initializer=proc_initializer, initargs=(cfg, _mapped_rxns), mp_context=context) as executor:
            print("Processing reactions w/ context: ", executor._mp_context)
            analyzed_reactions = list(
                tqdm(
                    executor.map(process_reaction, pred_rxns_to_do.iter_rows(named=True), chunksize=chunksize),
                    total=len(pred_rxns_to_do),
                    desc="Procesing reactions"
                )
            )
            
        analyzed_reactions = pl.from_dicts(analyzed_reactions, schema=predicted_reactions_schema)
        existing_reactions = pl.read_parquet("predicted_reactions.parquet")

        updated_reactions = update_table(existing_reactions, analyzed_reactions)
        updated_reactions.write_parquet("predicted_reactions.parquet")
        del analyzed_reactions, existing_reactions, updated_reactions, _mapped_rxns, pred_rxns_to_do

    # Retrieve info to update path stats
    path_ids_to_do = pl.scan_parquet("path_stats.parquet").filter(
        pl.col("mean_max_rxn_sim").is_null() | pl.col("mean_mean_rxn_sim").is_null() | pl.col("min_max_rxn_sim").is_null() | pl.col("min_mean_rxn_sim").is_null() | pl.col("feasibility_frac").is_null()
    ).select(
        pl.col("id")
    ).collect()["id"].to_list()

    paths_to_do = pl.scan_parquet("paths.parquet").filter(
        pl.col("path_id").is_in(path_ids_to_do)
    ).select(
        pl.col("path_id"),
        pl.col("rxn_id"),
        pl.col("rxn_type")
    ).collect()

    req_prids = paths_to_do.filter(
        pl.col("rxn_type") == "predicted"
    )["rxn_id"].explode().unique().to_list()

    req_rxns = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("id").is_in(req_prids)
    ).select(
        pl.col("id"),
        pl.col("dxgb_label"),
        pl.col("rxn_sims"),
    ).collect()

    # Known reactions get perfect scores
    req_krids = paths_to_do.filter(
        pl.col("rxn_type") == "known"
    )["rxn_id"].explode().unique().to_list()

    if len(req_krids) > 0:
        req_krs = pl.from_dicts(
            data=[{"id": rid, "dxgb_label": 1, "rxn_sims": [1]} for rid in req_krids],
            schema=req_rxns.schema
        )

        req_rxns = pl.concat((req_rxns, req_krs))
    
    paths_to_do = paths_to_do.group_by('path_id').agg(pl.col("rxn_id"))
    
    # Add reaction-derived summary stats to paths
    analyzed_path_stats = []
    for row in tqdm(paths_to_do.iter_rows(named=True), total=len(paths_to_do), desc="Updating path stats"):
        # If rxn_sims (& analogue_ids) for a predicted rxn is empty, there are 
        # no analogues. When calculating summary stats, treat this as rxn_sims = [0.0]
        # instead of skipping it for example.

        rxns = req_rxns.filter(
            pl.col("id").is_in(row['rxn_id'])
        ).with_columns(
            pl.col("rxn_sims").list.mean().fill_null(0).alias("mean_rxn_sims"),
            pl.col("rxn_sims").list.max().fill_null(0).alias("max_rxn_sims"),
        )

        analyzed_path_stats.append(
            {
                "id": row['path_id'],
                "feasibility_frac": rxns["dxgb_label"].mean(),
                "mean_max_rxn_sim": rxns["max_rxn_sims"].mean(),
                "mean_mean_rxn_sim": rxns["mean_rxn_sims"].mean(),
                "min_max_rxn_sim": rxns["max_rxn_sims"].min(),
                "min_mean_rxn_sim": rxns["mean_rxn_sims"].min(),
            }
        )

    analyzed_path_stats = pl.from_dicts(analyzed_path_stats, schema=path_stats_schema)
    existing_path_stats = pl.read_parquet("path_stats.parquet")
    updated_path_stats = update_table(existing_path_stats, analyzed_path_stats)
    updated_path_stats.write_parquet("path_stats.parquet")

if __name__ == "__main__":
    main()