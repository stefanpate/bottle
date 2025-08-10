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

rxn_type = pl.Enum(categories=["predicted", "known"])

paths_schema = pl.Schema(
    {
        "path_id": pl.String,
        "rxn_id": pl.String,
        "main_pdt_id": pl.String,
        "rxn_type": rxn_type,
        "generation": pl.Int32,
    }
)

path_stats_schema = pl.Schema(
    {
        "id": pl.String,
        "starters": pl.List(pl.String),
        "targets": pl.List(pl.String),
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