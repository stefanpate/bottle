import polars as pl

predicted_reactions_schema = pl.Schema(
    {
        "id": pl.String,
        "smarts": pl.String,
        "am_smarts": pl.String,
        "dxgb_label": pl.Int32,
        "rxn_sims": pl.List(pl.Float32),
        "analogue_ids": pl.List(pl.String),
        "rules": pl.List(pl.String),
    }
)

found_paths_schema = pl.Schema(
    {
        "id": pl.String,
        "starters": pl.List(pl.String),
        "targets": pl.List(pl.String),
        "reactions": pl.List(pl.String),
        "dg_opt": pl.List(pl.Float32),
        "dg_err": pl.List(pl.Float32),
        "starter_ids": pl.List(pl.String),
        "target_ids": pl.List(pl.String),
        "mdf": pl.Float32,
        "mean_max_rxn_sim": pl.Float32,
        "mean_mean_rxn_sim": pl.Float32,
        "min_max_rxn_sim": pl.Float32,
        "min_mean_rxn_sim": pl.Float32,
        "feasibility_frac": pl.Float32,
    }
)

expansion_reactions_schema = pl.Schema(
    {
        "am_smarts": pl.String,
        "rules": pl.List(pl.String),
    }
)