import hydra
from omegaconf import DictConfig
import polars as pl
from pathlib import Path
from src.chem_draw import draw_reaction
from time import perf_counter

@hydra.main(version_base=None, config_path="../conf", config_name="draw_reactions")
def main(cfg: DictConfig) -> None:

    if not Path("paths.parquet").exists():
        return
    
    if not Path("svgs").exists():
        Path("svgs").mkdir()
    
    existing_drawn_reactions = [kr.stem for kr in Path("svgs").glob("*.svg")]
    
    rids_to_draw = pl.scan_parquet("paths.parquet").unique(subset=["rxn_id"]).filter(
        ~pl.col("rxn_id").is_in(existing_drawn_reactions)
    ).select(
        pl.col("rxn_id")
    ).collect()['rxn_id'].to_list()

    pr_smarts = pl.scan_parquet("predicted_reactions.parquet").filter(
        pl.col("id").is_in(rids_to_draw)
    ).select(
        pl.col("id"),
        pl.col("smarts"),
    ).collect()

    kr_smarts = pl.scan_parquet(cfg.filepaths.known_reactions).filter(
        pl.col("id").is_in(rids_to_draw)
    ).select(
        pl.col("id"),
        pl.col("smarts"),
    ).collect()

    rxns_to_draw = pl.concat([pr_smarts, kr_smarts]).unique(subset=["id"])

    tic = perf_counter()
    print(f"Drawing {len(rxns_to_draw)} known reactions")
    for row in rxns_to_draw.iter_rows(named=True):
        rxn = draw_reaction(row['smarts'], auto_scl=True)
        rxn.save(f"svgs/{row['id']}.svg")
    toc = perf_counter()
    print(f"Reaction drawing completed in {toc - tic : .2f} seconds")

if __name__ == "__main__":
    main()