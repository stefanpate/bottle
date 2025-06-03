import hydra
from omegaconf import DictConfig, OmegaConf
from tqdm import tqdm
import polars as pl
from concurrent.futures import ProcessPoolExecutor
from hydra.utils import instantiate
from functools import partial
from ergochemics.mapping import get_reaction_center
from pathlib import Path
from minedatabase.pickaxe import Pickaxe
from rdkit import Chem
from src.post_processing import Expansion
import numpy as np
from time import perf_counter

def dxgb_initializer(cfg: DictConfig):
    global dxgb, mfper, _fingerprint
    dxgb = instantiate(cfg.dxgb)
    mfper = instantiate(cfg.mfper)
    _fingerprint = partial(mfper.fingerprint, rc_dist_ub=cfg.rc_dist_ub)

def get_lhs_block_rc(am_smarts: str) -> list[int]:
    return get_reaction_center(am_smarts, mode="combined")[0]

def extract_rule_names(fwd_expansions: list[str], reverse_expansions: list[str]) -> str:
    rule_names = []
    for fwd, rev in zip(fwd_expansions, reverse_expansions):
        if fwd is not None and rev is not None:
            fwd_rules = fwd.split("_rules_")[1]
            rev_rules = rev.split("_rules_")[1]

            if fwd_rules != rev_rules:
                raise ValueError(f"Forward and reverse expansions have different rules: {fwd_rules} != {rev_rules}")
            rules = fwd_rules
        elif rev is None:
            rules = fwd.split("_rules_")[1]
        elif fwd is None:
            rules = rev.split("_rules_")[1]
        else:
            raise ValueError("Both forward and reverse expansions are None")
        rule_names.append(f"{rules}_rules")
    return rule_names

OmegaConf.register_new_resolver("extract_rule_names", extract_rule_names)

# TODO: modify for this use case
def process_reactions(pk: Pickaxe, mapped_rxn: pl.DataFrame) -> float:
    # Collect starters
    starters = {}
    for v in pk.compounds.values():
        if v["Type"].startswith("Start"):
            starters[v['_id']] = v["ID"]

    # Calculate morgan fps for all mapped reactions
    mapped_rxn["mol"] = mapped_rxn["smarts"].apply(lambda x : Chem.MolFromSmiles(x.split('>>')[0]))
    mapped_rxn["mfp"] = mapped_rxn.apply(lambda x : _fingerprint(x.mol, x.reaction_center), axis=1)

    data = []
    for v in pk.reactions.values():

        query_smarts = v["Operator_aligned_smarts"]
        query_am_smarts = v["am_rxn"]
        query_lhs_mol = Chem.MolFromSmiles(query_smarts.split('>>')[0])
        query_lhs_block_rc = get_lhs_block_rc(query_am_smarts)

        if query_lhs_mol is None:
            continue

        rules = set([int(elt.split('_')[0]) for elt in v["Operators"]])
        analogues = mapped_rxn.loc[mapped_rxn.rule_id.isin(rules)]

        if analogues.empty:
            max_sim = 0.0
            nearest_kr = ''
            nearest_krid = ''
        else:
            mfps = np.vstack(analogues["mfp"])
            query_mfp = _fingerprint(
                mol=query_lhs_mol,
                reaction_center=query_lhs_block_rc
            ).reshape(-1, 1)
            sims = np.matmul(mfps, query_mfp) / (np.linalg.norm(mfps, axis=1).reshape(-1, 1) * np.linalg.norm(query_mfp))
            sims = sims.reshape(-1,)
            max_sim = float(np.max(sims))
            max_idx = int(np.argmax(sims))
            nearest_kr = analogues.iloc[max_idx].smarts
            nearest_krid = analogues.iloc[max_idx].rxn_id

        is_feasible = dxgb.predict_label(query_smarts)

        data.append(
            [
                v['_id'],
                query_smarts,
                query_am_smarts,
                is_feasible,
                max_sim,
                nearest_kr,
                nearest_krid,
                list(v['Operators'])
            ]
        )

    return data

# TODO: write custom resolver for mapped reaction paths

@hydra.main(version_base=None, config_path="../conf", config_name="mechanism_processing")
def main(cfg: DictConfig) -> None:

    expansions = []
    for fwd, rev in zip(cfg.forward_expansions, cfg.reverse_expansions):
        pk = Expansion(
            forward=Path(cfg.filepaths.raw_expansions) / fwd if fwd else None,
            reverse=Path(cfg.filepaths.raw_expansions) / rev if rev else None,
        )
        print("Searching for paths")
        tic = perf_counter()
        paths = pk.find_paths()
        toc = perf_counter()
        print(f"Found {sum([len(v) for v in paths.values()])} paths in  {toc - tic : .2f} seconds")
        pk.prune(paths)
        expansions.append(pk)
        print(f"Pruned expansion to {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
        
    # Required to look up analogues of reveresed reactions from retro expansions
    rxn_reverses = pl.scan_parquet(
        Path(cfg.filepaths.known_reactions)
    ).select(
        pl.col("id"),
        pl.col("reverse")
    ).collect()
    rxn_reverses = dict(rxn_reverses.iter_rows())
    
    mapped_rxns = []
    for rn in cfg.rule_names:
        df = pl.read_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / f"{cfg.krs_name}_x_{rn}.parquet"
        )
        df = df.with_columns(
            reaction_center=pl.col("am_smarts").map_elements(
                get_lhs_block_rc, return_dtype=pl.List(pl.Int64)
            )
        )
        mapped_rxns.append(df)

    # Process reactions
    with ProcessPoolExecutor(max_workers=len(expansions), initializer=dxgb_initializer, initargs=(cfg,)) as executor:
        results = list(
            tqdm(
                executor.map(process_reactions, expansions, mapped_rxns, chunksize=1),
                total=len(expansions),
                desc="Procesing reactions"
            )
        )

    # Save reaction metrics
    columns = ["id", "smarts", "am_smarts", "dxgb_label", "max_rxn_sim", "nearest_analogue", "nearest_analogue_id", "rules"] 
    dfes = []
    for exp, res in zip(cfg.expansion_fns, results):
        df = pd.DataFrame(data=res, columns=columns)
        df["expansion"] = exp
        dfes.append(df)

    full_df = pd.concat(dfes)
    full_df.to_parquet(f"{cfg.expansion_name}_reaction_metrics.parquet")

if __name__ == "__main__":
    main()