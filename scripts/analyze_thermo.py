import hydra
from omegaconf import DictConfig
from tqdm import tqdm
import polars as pl
from pathlib import Path
from ergochemics.standardize import hash_molecule
from time import perf_counter
from logging import getLogger
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_api.phased_reaction import PhasedReaction
from equilibrator_api import Q_, ComponentContribution
import equilibrator_assets.chemaxon as _chemaxon
import cvxpy

from src.schemas import path_stats_schema
from src.post_processing import pick_constraints_for_MDF
from src.pka_plugins import MolGPKA

# Monkey patch pka predictor
_predictor = MolGPKA()
_chemaxon.get_dissociation_constants = _predictor.get_dissociation_constants

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
    ).select(
        pl.col('id')
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
            "coco_id": [hash_molecule(smi) for smi in unique_cpds],
            "name": [cid2name.get(smi, "unknown") for smi in unique_cpds],
        }
    ).to_pandas()

    logger.info(f"{len(unique_cpds)} unique compounds / {len(rxns_to_do)} reactions / {len(path_stats_to_do)} paths to analyze")

    logger.info("Adding compounds to local cache...")
    lc = LocalCompoundCache(ccache_path=cfg.eq_cache)
    start = perf_counter()
    lc.add_compounds(
        unique_cpds,
        mol_format="smiles",
        bypass_chemaxon=True, # Still adds compound even if cxcalc fails
        save_empty_compounds=True,
    )
    end = perf_counter()
    logger.info(f"Added {len(unique_cpds)} compounds to local cache in {end - start:.2f} seconds")

    # Gather required compounds
    eq_cpds = lc.get_compounds(unique_cpds['struct'].to_list())
    smiles2cid = dict(zip(unique_cpds['struct'].to_list(), unique_cpds['coco_id'].to_list()))
    cid_2_eq_cpd = dict(zip(unique_cpds['coco_id'].to_list(), [elt.compound for elt in eq_cpds]))

    # Gather required reactions
    eq_cpd_getter = lambda cid: cid_2_eq_cpd.get(cid, None)
    rid_2_eq_rxn = {}
    for row in rxns_to_do.iter_rows(named=True):
        lhs, rhs = [elt.split('.') for elt in row['smarts'].split(">>")]
        lhs = [smiles2cid[elt] for elt in lhs]
        rhs = [smiles2cid[elt] for elt in rhs]
        lhs = " + ".join(lhs)
        rhs = " + ".join(rhs)
        rxn_string = f"{lhs} = {rhs}"
        rid_2_eq_rxn[row['id']] = PhasedReaction.parse_formula(eq_cpd_getter, rxn_string)

    # Calculate path mdfs
    cc = ComponentContribution(ccache=lc.ccache)
    cc.p_h = Q_(cfg.p_h)
    cc.p_mg = Q_(cfg.p_mg)
    cc.ionic_strength = Q_(cfg.ionic_strength)
    cc.temperature = Q_(cfg.temperature)

    analyzed_path_stats = []
    for row in tqdm(path_stats_to_do.iter_rows(named=True), total=len(path_stats_to_do), desc="Calculating path MDFs"):
        rids = paths_to_do.filter(
            pl.col("path_id") == row['id']
        ).sort(
            pl.col("generation"),
            descending=False
        )['rxn_id'].to_list()
        
        failed_retrieval = False
        path_eq_rxns = []
        for rid in rids:
            eq_rxn = rid_2_eq_rxn.get(rid, None)
            
            if eq_rxn is None:
                failed_retrieval = True
                break
            
            path_eq_rxns.append(eq_rxn)
        
        if failed_retrieval:
            logger.warning(f"Failed to retrieve all reactions for path {row['id']}, skipping MDF calculation")
            continue

        standard_dgr_prime, standard_dgr_uncertainty = cc.standard_dg_prime_multi(
            path_eq_rxns,
            uncertainty_representation="fullrank"
        )

        S = cc.create_stoichiometric_matrix_from_reaction_objects(path_eq_rxns)
        Nc, Nr = S.shape
        RT = cc.RT

        ln_conc = cvxpy.Variable(
            shape=Nc, name="metabolite log concentration"
        )
        B = cvxpy.Variable()
        dg_prime = -(
            standard_dgr_prime.m_as("kJ/mol") + RT.m_as("kJ/mol") * S.values.T @ ln_conc
        )
        
        constraints = pick_constraints_for_MDF(S, Nc, Nr, ln_conc, dg_prime, B, cfg.conc_lb, cfg.conc_ub)

        # Solve the MDF problem
        prob_max = cvxpy.Problem(cvxpy.Maximize(B), constraints)
        prob_max.solve()

        if prob_max.value is None:
            logger.warning(f"Path {row['id']} MDF optimization failed with status {prob_max.status}, skipping")
            continue

        # Convert the results to the correct units
        dg_opt = [
            Q_(val, "kilojoule/mole").magnitude for val in list(-dg_prime.value)
        ]

        dg_err = [
            unc[0].magnitude for unc in standard_dgr_uncertainty
        ]

        analyzed_path_stats.append(
            {
                "id": row['id'],
                "mdf": prob_max.value,
                "dg_opt": dg_opt,
                "dg_err": dg_err,
            }
        )

    analyzed_path_stats = pl.from_dicts(analyzed_path_stats, schema=path_stats_schema)
    existing_path_stats = pl.read_parquet("path_stats.parquet")
    updated_path_stats = update_table(existing_path_stats, analyzed_path_stats)
    updated_path_stats.write_parquet("path_stats.parquet")

if __name__ == "__main__":
    main()