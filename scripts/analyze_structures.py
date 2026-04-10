import hydra
from omegaconf import DictConfig
from tqdm import tqdm
import polars as pl
from concurrent.futures import ProcessPoolExecutor
from hydra.utils import instantiate
from ergochemics.mapping import get_reaction_center
from ergochemics.similarity import ReactionFingerprinter, MolFeaturizer
from pathlib import Path
from rdkit import Chem
from src.schemas import path_stats_schema, predicted_reactions_schema
import numpy as np
from typing import Any
from multiprocessing import get_context
from logging import getLogger
from time import perf_counter

def update_table(existing: pl.DataFrame, analyzed: pl.DataFrame, on: str = "id") -> pl.DataFrame:
    ''' Update existing table with analyzed data, filling in missing values.'''
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

def proc_initializer(dxgb_cfg: DictConfig, _am_krs: pl.DataFrame, _S: np.ndarray) -> None:    
    print("Initializing process for reaction mapping", flush=True)
    global dxgb, am_krs, S
    dxgb = instantiate(dxgb_cfg)
    am_krs = _am_krs
    S = _S 

def process_reaction(reaction: dict[str, Any]) -> dict[str, Any]:
    ''' Process a single reaction: compute similarity scores to known reactions and predict feasibility.'''
    query_smarts = reaction["smarts"]
    query_index = reaction["index"]
    
    sims = S[query_index, :]
    srt_sims = sorted(sims, reverse=True)
    srt_idxs = np.argsort(sims)[::-1]
    first_zero_idx = np.argmax(np.array(srt_sims) == 0.0)

    if first_zero_idx == 0:
        srt_sims = []
        srt_krids = []
    else:
        srt_sims = srt_sims[:first_zero_idx]
        srt_idxs = srt_idxs[:first_zero_idx]
        srt_krids = am_krs['rxn_id'][srt_idxs].to_list()

        # De-duplicate, keeping the highest similarity score for each known reaction
        seen_krids = set()
        unique_srt_sims = []
        unique_srt_krids = []
        for sim, krid in zip(srt_sims, srt_krids):
            if krid not in seen_krids:
                unique_srt_sims.append(sim)
                unique_srt_krids.append(krid)
                seen_krids.add(krid)
        srt_sims = unique_srt_sims
        srt_krids = unique_srt_krids

    is_feasible = dxgb.predict_label(query_smarts)
    reaction["dxgb_label"] = is_feasible
    reaction["rxn_sims"] = srt_sims
    reaction["analogue_ids"] = srt_krids
    return reaction

def get_rc_patts(am_rxn: str) -> tuple[str, str]:
    ''' Get reaction center pattern SMILES (left and right) from an atom-mapped reaction.'''
    sort_side = lambda side: ".".join(sorted(side.split(".")))
    lrc_idx, rrc_idx = get_reaction_center(am_rxn, mode="combined")
    lhs, rhs = [Chem.MolFromSmiles(smi) for smi in am_rxn.split(">>")]

    for mol in [lhs, rhs]:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

    lrc = Chem.MolFragmentToSmiles(lhs, lrc_idx)
    rrc = Chem.MolFragmentToSmiles(rhs, rrc_idx)

    lrc = sort_side(lrc)
    rrc = sort_side(rrc)
    return (lrc, rrc)

def tani_sim(X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
    ''' Compute Tanimoto similarity between two sets of fingerprints.'''
    A = X1 @ X2.T
    norms1 = X1.sum(axis=1).reshape(-1, 1)
    norms2 = X2.sum(axis=1).reshape(1, -1)
    N = norms1 + norms2
    S = A / (N - A)
    return S

def rcmfp_S(left: pl.DataFrame, right: pl.DataFrame, use_rc: bool, filter_by_rc: bool, radius: int, length: int, rc_dist_ub: int) -> np.ndarray:
    ''' 
    Compute RCMFP similarity matrix between the peirs of reactions 
    from the left and right dataframes.
    
    Args
    ----
    left: pd.DataFrame
        DataFrame containing the reactions with 'am_smarts' column
    right: pd.DataFrame
        DataFrame containing the reactions with 'am_smarts' column
    use_rc: bool
        Whether to use reaction center information in fingerprinting.
    filter_by_rc: bool
        Whether to filter similarity scores by reaction center patterns.
    radius: int
        Radius for fingerprinting.
    length: int
        Length of the fingerprint.
    rc_dist_ub: int
        Upper bound on bondwise distance from reaction center to consider in fingerprinting.

    Returns
    -------
    S: np.ndarray
        Similarity matrix of dimensions (len(left) x len(right))
    '''
    rxn_fper = ReactionFingerprinter(
        radius=radius,
        length=length,
        mol_featurizer=MolFeaturizer(),
    )

    left_rc_patts = []
    left_fps = []
    i = 0
    for row in left.iter_rows(named=True):
        smi = row['am_smarts']
        rc_patt = get_rc_patts(smi)
        left_rc_patts.append(rc_patt)
        fp = rxn_fper.fingerprint(smi, use_rc=use_rc, rc_dist_ub=rc_dist_ub)
        left_fps.append(fp)
        i += 1

    left_fps = np.vstack(left_fps).astype(np.float32)
    
    right_rc_patts = []
    right_fps = []
    i = 0
    for row in right.iter_rows(named=True):
        smi = row['am_smarts']
        rc_patt = get_rc_patts(smi)
        right_rc_patts.append(rc_patt)
        fp = rxn_fper.fingerprint(smi, use_rc=use_rc, rc_dist_ub=rc_dist_ub)
        right_fps.append(fp)
        i += 1

    right_fps =  np.vstack(right_fps).astype(np.float32)
    d = right_fps.shape[1]
    right_fps_rev = np.concatenate([right_fps[:, d//2:], right_fps[:, :d//2]], axis=1)

    S1 = tani_sim(left_fps, right_fps)
    S2 = tani_sim(left_fps, right_fps_rev)
    S = np.maximum(S1, S2)

    if filter_by_rc:
        filtered_S = np.zeros_like(S)
        for i, patt_i in enumerate(left_rc_patts):
            for j, patt_j in enumerate(right_rc_patts):
                if patt_i == patt_j:
                    filtered_S[i, j] = S[i, j]
                    continue

                if tuple(reversed(patt_i)) == patt_j:
                    filtered_S[i, j] = S[i, j]
                    continue

        return filtered_S
    else:
        return S

logger = getLogger(__name__)

@hydra.main(version_base=None, config_path="../conf", config_name="analyze_structures")
def main(cfg: DictConfig) -> None:

    if not Path("paths.parquet").exists():
        return
    
    logger.info("Loading predicted reactions to analyze...")

    pred_rxns_to_do = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("dxgb_label").is_null() | pl.col("rxn_sims").is_null() | pl.col("analogue_ids").is_null()
    ).collect().with_row_index()

    if not len(pred_rxns_to_do) == 0:
        logger.info(f"Found {len(pred_rxns_to_do)} predicted reactions to analyze.")
        logger.info("Loading atom-mapped known reactions...")

        mechinformed_am_krs = pl.scan_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / "mapped_known_reactions_x_mechinformed_rules.parquet",
        ).select(
            pl.col("rxn_id"),
            pl.col("am_smarts")
        ).collect()
        
        rc_0_am_krs = pl.scan_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / "mapped_known_reactions_x_rc_plus_0_rules.parquet",
        ).select(
            pl.col("rxn_id"),
            pl.col("am_smarts")
        ).collect()
        
        _am_krs = pl.concat([rc_0_am_krs, mechinformed_am_krs]).with_row_index()

        logger.info(f"Computing reaction similarities between {len(pred_rxns_to_do)} predicted reactions and {len(_am_krs)} known reactions.")
        tic = perf_counter()
        _S = rcmfp_S(
            left=pred_rxns_to_do,
            right=_am_krs,
            use_rc=cfg.use_rc_in_fingerprints,
            filter_by_rc=cfg.filter_by_rc,
            radius=cfg.mfper.radius,
            length=cfg.mfper.length,
            rc_dist_ub=cfg.rc_dist_ub
        )
        toc = perf_counter()
        logger.info(f"Computed similarity matrix in {toc - tic:.2f} seconds.")
        
        # Process reactions
        chunksize = max(1, int(len(pred_rxns_to_do) / cfg.processes))
        context = get_context("spawn") # Polars hates fork
        logger.info(f"Processing reactions using {cfg.processes} processes with chunksize {chunksize}...")
        with ProcessPoolExecutor(max_workers=cfg.processes, initializer=proc_initializer, initargs=(cfg.dxgb, _am_krs, _S), mp_context=context) as executor:
            print("Processing reactions w/ context: ", executor._mp_context)
            analyzed_reactions = list(
                tqdm(
                    executor.map(process_reaction, pred_rxns_to_do.iter_rows(named=True), chunksize=chunksize),
                    total=len(pred_rxns_to_do),
                    desc="Procesing reactions"
                )
            )
            
        logger.info("Updating predicted reactions table...")
        analyzed_reactions = pl.from_dicts(analyzed_reactions, schema=predicted_reactions_schema)
        existing_reactions = pl.read_parquet("predicted_reactions.parquet")

        updated_reactions = update_table(existing_reactions, analyzed_reactions)
        updated_reactions.write_parquet("predicted_reactions.parquet")
        
        del analyzed_reactions, existing_reactions, updated_reactions, _am_krs, pred_rxns_to_do, _S

    # Retrieve info to update path stats
    logger.info("Updating path statistics...")
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