import hydra
from omegaconf import DictConfig, OmegaConf
from tqdm import tqdm
import polars as pl
from concurrent.futures import ProcessPoolExecutor
from hydra.utils import instantiate
from functools import partial
from itertools import chain
from ergochemics.mapping import get_reaction_center
from pathlib import Path
from rdkit import Chem
from src.post_processing import Expansion, get_path_id
from src.schemas import found_paths_schema, predicted_reactions_schema
from src.chem_draw import draw_reaction
import numpy as np
from time import perf_counter
from typing import Any
from multiprocessing import get_context
import json

def proc_initializer(cfg: DictConfig) -> None:
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

    # Create helper dfs
    mapped_rxns = []
    unique_rule_names = set(chain(*cfg.rule_names))
    unique_rule_names.discard(None)  # Remove None if present
    for rn in unique_rule_names:
        df = pl.read_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / f"mapped_known_reactions_x_{rn}.parquet"
        )
        df = df.with_columns(
            reaction_center=pl.col("am_smarts").map_elements(
                get_lhs_block_rc, return_dtype=pl.List(pl.Int64)
            ),
            rule_id=pl.col("rule_id").map_elements(
                lambda x: f"{rn}:{x}",
                return_dtype=pl.String
            )
        )
        df = df.with_columns(
            pl.struct(["smarts", "reaction_center"]).map_elements(
                lambda x: sma_rc_fp(x["smarts"], x["reaction_center"]),
                return_dtype=pl.Object
            ).alias("mfp")
        )
        mapped_rxns.append(df)

    mapped_rxns = pl.concat(mapped_rxns)
    
    # Required to look up analogues of reveresed reactions from retro expansions
    rxn_reverses = pl.scan_parquet(
        Path(cfg.filepaths.known_reactions)
    ).select(
        pl.col("id"),
        pl.col("reverse")
    ).collect()
    rxn_reverses = dict(rxn_reverses.iter_rows())

def get_lhs_block_rc(am_smarts: str) -> list[int]:
    return get_reaction_center(am_smarts, mode="combined")[0]

def extract_rule_names(fwd_expansions: list[str], reverse_expansions: list[str]) -> str:
    rule_names = []
    for fwd, rev in zip(fwd_expansions, reverse_expansions):
        if fwd is not None and rev is not None:
            fwd_rules = fwd.split("_rules_")[1]
            rev_rules = rev.split("_rules_")[1]
            rule_names.append((f"{fwd_rules}_rules", f"{rev_rules}_rules"))
        elif rev is None:
            rules = fwd.split("_rules_")[1]
            rule_names.append((f"{rules}_rules", None))
        elif fwd is None:
            rules = rev.split("_rules_")[1]
            rule_names.append((None, f"{rules}_rules"))
        else:
            raise ValueError("Both forward and reverse expansions are None")
 
    return rule_names

OmegaConf.register_new_resolver("extract_rule_names", extract_rule_names)

def process_reaction(reaction: dict[str, Any]) -> dict[str, Any]:
    reaction = {k: v for k, v in reaction.items()} # Defensive copy
    query_smarts = reaction["smarts"]
    query_am_smarts = reaction["am_smarts"]
    query_lhs_mol = Chem.MolFromSmiles(query_smarts.split('>>')[0])
    query_lhs_block_rc = get_lhs_block_rc(query_am_smarts)

    if query_lhs_mol is None:
        return reaction
    
    # In case of retro expansion, must get reactions of the reverse direction
    if reaction["reversed"]:
        rev_krs = mapped_rxns.filter(
            pl.col("rule_id").is_in(reaction["rules"])
        )["rxn_id"]
        fwd_krs = [rxn_reverses[kr] for kr in rev_krs if kr in rxn_reverses]
        analogues = mapped_rxns.filter(
            pl.col("rxn_id").is_in(fwd_krs)
        )
    else: # For forward expansions, we can use the rules to look up analogues directly
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
    reaction.pop("reversed", None) # Don't want to save this

    return reaction

@hydra.main(version_base=None, config_path="../conf", config_name="mechanism_processing")
def main(cfg: DictConfig) -> None:
    # Check for existing paths and reactions
    if Path("predicted_reactions.parquet").exists():
        existing_reactions = pl.scan_parquet(
            "predicted_reactions.parquet"
        ).select(
            pl.col("id")
        ).collect()
    else:
        existing_reactions = pl.DataFrame(schema=predicted_reactions_schema)
        existing_reactions.write_parquet("predicted_reactions.parquet")

    if Path("found_paths.parquet").exists():
        existing_paths = pl.scan_parquet(
            "found_paths.parquet"
        ).select(
            pl.col("id")
        ).collect()
    else:
        existing_paths = pl.DataFrame(schema=found_paths_schema)
        existing_paths.write_parquet("found_paths.parquet")
    
    predicted_reactions = {}
    paths = {}
    for fwd, rev, (fwd_rules_name, rev_rules_name) in zip(cfg.forward_expansions, cfg.reverse_expansions, cfg.rule_names):

        pk = Expansion(
            forward=Path(cfg.filepaths.raw_expansions) / fwd if fwd else None,
            reverse=Path(cfg.filepaths.raw_expansions) / rev if rev else None,
        )
        
        # Find paths
        print("Searching for paths")
        tic = perf_counter()
        this_paths = pk.find_paths()
        toc = perf_counter()
        print(f"Path finding completed in  {toc - tic : .2f} seconds")
        tic = perf_counter()
        print("Pruning paths")
        pk.prune(this_paths)
        toc = perf_counter()
        
        print(f"Path pruning completed in {toc - tic : .2f} seconds")
        if not this_paths:
            print(f"No paths found for {fwd} and {rev}. Skipping.")
            continue
        else:
            print(f"Found {len(this_paths)} paths for {fwd} and {rev}")
        
        for sids, tids, rids in this_paths:
            path_id = get_path_id(rids)

            if path_id in existing_paths['id']:
                continue

            # Entering some default values here, to be updated later, 
            # in order to adhere to the schema
            paths[path_id] = {
                "id": path_id,
                "starters": [pk.starters[sid] for sid in sids],
                "targets": [pk.targets[tid] for tid in tids],
                "reactions": rids,
                "dg_opt": None,
                "dg_err": None,
                "starter_ids": sids,
                "target_ids": tids,
                "mdf": None,
                "mean_max_rxn_sim": 0.0,
                "mean_mean_rxn_sim": 0.0,
                "min_max_rxn_sim": 0.0,
                "min_mean_rxn_sim": 0.0,
                "feasibility_frac": 0.0,
            }
        
        # Collect predicted reactions post-pruning
        for k, v in pk.reactions.items():

            if k in existing_reactions["id"]:
                continue

            is_reversed = v["reversed"]
            rule_name = rev_rules_name if is_reversed else fwd_rules_name
            predicted_reactions[k] = {
                "id": k,
                "smarts": v["Operator_aligned_smarts"],
                "am_smarts": v["am_rxn"],
                "dxgb_label": None,
                "rxn_sims": None,
                "analogue_ids": None,
                "rules": [f"{rule_name}:{elt.split('_')[0]}" for elt in v["Operators"]],
                "reversed": is_reversed,
            }

    # Process reactions
    predicted_reactions = list(predicted_reactions.values())
    chunksize = max(1, int(len(predicted_reactions) / cfg.processes))
    context = get_context("spawn") # Polars hates fork
    with ProcessPoolExecutor(max_workers=cfg.processes, initializer=proc_initializer, initargs=(cfg,), mp_context=context) as executor:
        print("Processing reactions w/ context: ", executor._mp_context)
        analyzed_reactions = list(
            tqdm(
                executor.map(process_reaction, predicted_reactions, chunksize=chunksize),
                total=len(predicted_reactions),
                desc="Procesing reactions"
            )
        )
        
    analyzed_reactions = pl.from_dicts(analyzed_reactions, schema=predicted_reactions_schema)

    # Add reaction-derived summary stats to paths
    for id, path in paths.items():           
        rxn = analyzed_reactions.filter(pl.col("id").is_in(path["reactions"])).select(
            pl.col("dxgb_label"),
            pl.col("rxn_sims"),
        )

        # If rxn_sims (& analogue_ids) for a predicted rxn is empty, there are 
        # no analogues. When calculating summary stats, treat this as rxn_sims = [0.0]
        # instead of skipping it for example.
        paths[id]["feasibility_frac"] = rxn["dxgb_label"].mean()
        paths[id]["mean_max_rxn_sim"] = np.mean([sims.max() or 0 for sims in rxn["rxn_sims"]] or [0])
        paths[id]["mean_mean_rxn_sim"] = np.mean([sims.mean() or 0 for sims in rxn["rxn_sims"]] or [0])
        paths[id]["min_max_rxn_sim"] = np.min([sims.max() or 0 for sims in rxn["rxn_sims"]] or [0])
        paths[id]["min_mean_rxn_sim"] = np.min([sims.mean() or 0 for sims in rxn["rxn_sims"]] or [0])
    
    # Save paths
    paths = pl.from_dicts(list(paths.values()), schema=found_paths_schema)
    existing_paths = pl.read_parquet("found_paths.parquet")
    paths = pl.concat((existing_paths, paths))
    paths.write_parquet("found_paths.parquet")

    # Generate reaction images
    if not Path("svgs").exists():
        Path("svgs").mkdir()
    
    if not Path("svgs/known").exists():
        Path("svgs/known").mkdir()
        existing_kr_svgs = set()
    else:
        existing_kr_svgs = set([fn.name.removesuffix(".svg") for fn in Path("svgs/known").glob("*.svg")])
    
    if not Path("svgs/predicted").exists():
        Path("svgs/predicted").mkdir()
            
    krids_to_draw = set(chain(*analyzed_reactions["analogue_ids"])) - existing_kr_svgs
    krs_to_draw = pl.scan_parquet(
        Path(cfg.filepaths.known_reactions)
    ).filter(
        pl.col("id").is_in(krids_to_draw)
    ).select(
        pl.col("id"),
        pl.col("smarts")
    ).collect()


    tic = perf_counter()
    print(f"Drawing {len(krs_to_draw)} known reactions")
    for row in krs_to_draw.iter_rows(named=True):
        rxn = draw_reaction(row['smarts'], auto_scl=True)
        rxn.save(f"svgs/known/{row['id']}.svg")
    toc = perf_counter()
    print(f"Drawing known reactions completed in {toc - tic : .2f} seconds")

    tic = perf_counter()
    print(f"Drawing {len(analyzed_reactions)} predicted reactions")
    for row in analyzed_reactions.iter_rows(named=True):
        rxn = draw_reaction(row['smarts'], auto_scl=True)
        rxn.save(f"svgs/predicted/{row['id']}.svg")
    toc = perf_counter()
    print(f"Drawing predicted reactions completed in {toc - tic : .2f} seconds")
    
    # Save predicted reactions
    existing_reactions = pl.read_parquet("predicted_reactions.parquet")
    analyzed_reactions = pl.concat((analyzed_reactions, existing_reactions))
    analyzed_reactions.write_parquet("predicted_reactions.parquet")

if __name__ == "__main__":
    main()