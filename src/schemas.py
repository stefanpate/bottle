import polars as pl

predicted_reactions_schema = {
    "id": pl.String,
    "smarts": pl.String,
    "am_smarts": pl.String,
    "dxgb_label": pl.Int32,
    "rxn_sims": pl.List(pl.Float32),
    "analogue_ids": pl.List(pl.String),
    "rules": pl.List(pl.String), 
}

found_paths_schema = {
    "id": pl.String,
    "starter": pl.String,
    "target": pl.String,
    "reactions": pl.List(pl.String),
    "dg_opt": pl.List(pl.Float32),
    "dg_err": pl.List(pl.Float32),
    "starter_id": pl.String,
    "target_id": pl.String,
    "mdf": pl.Float32
}