import hydra
from omegaconf import DictConfig
import polars as pl
from pathlib import Path
from src.chem_draw import draw_reaction
from tqdm import tqdm

@hydra.main(version_base=None, config_path="../conf", config_name="draw_reactions")
def main(cfg: DictConfig) -> None:

    if not Path("paths.parquet").exists():
        return
    
    if not Path("svgs").exists():
        Path("svgs").mkdir()
    
    existing_drawn_reactions = [kr.stem for kr in Path("svgs").glob("*.svg")]
    
    path_rids = pl.scan_parquet("paths.parquet").unique(subset=["rxn_id"]).filter(
        ~pl.col("rxn_id").is_in(existing_drawn_reactions)
    ).select(
        pl.col("rxn_id")
    ).collect()['rxn_id'].to_list()

    prs = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("id").is_in(path_rids)
    ).select(
        pl.col("id"),
        pl.col("smarts"),
        pl.col("analogue_ids"),
    ).collect()

    pr_smarts = prs.select(
        pl.col("id"),
        pl.col("smarts"),
    )

    analogue_ids = [id for id in prs['analogue_ids'].explode().unique() if id not in existing_drawn_reactions]

    kr_smarts = pl.scan_parquet(cfg.filepaths.known_reactions).filter(
        pl.col("id").is_in(path_rids + analogue_ids)
    ).select(
        pl.col("id"),
        pl.col("smarts"),
    ).collect()

    rxns_to_draw = pl.concat([pr_smarts, kr_smarts]).unique(subset=["id"])

    for row in tqdm((rxns_to_draw.iter_rows(named=True)), total=len(rxns_to_draw), desc="Drawing reactions"):
        rxn = draw_reaction(row['smarts'], auto_scl=True)
        rxn.save(f"svgs/{row['id']}.svg")

if __name__ == "__main__":
    main()