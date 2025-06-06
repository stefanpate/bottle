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
        "reversed": pl.Boolean,
    }
)

found_paths_schema = pl.Schema(
    {
        "id": pl.String,
        "starter": pl.String,
        "target": pl.String,
        "reactions": pl.List(pl.String),
        "dg_opt": pl.List(pl.Float32),
        "dg_err": pl.List(pl.Float32),
        "starter_id": pl.String,
        "target_id": pl.String,
        "mdf": pl.Float32,
        "mean_max_rxn_sim": pl.Float32,
        "mean_mean_rxn_sim": pl.Float32,
        "min_max_rxn_sim": pl.Float32,
        "min_mean_rxn_sim": pl.Float32,
        "feasibility_frac": pl.Float32,
    }
)

known_compounds_schema = pl.Schema(
    {
        "id": pl.String,
        "smiles": pl.String,
        "names": pl.List(pl.String),
        "n_atoms": pl.Int32,
    }
)

existence_enum = pl.Enum(
    [
        "Predicted",
        "Uncertain",
        "Inferred from homology",
        "Evidence at transcript level",
        "Evidence at protein level",
    ]
)

reviewed_enum = pl.Enum(
    [
        "reviewed",
    ]
)

enzymes_schema = pl.Schema(
    {
        "id": pl.String,
        "sequence": pl.String,
        "existence": existence_enum,
        "reviewed": reviewed_enum,
        "ec": pl.String,
        "organism": pl.String,
        "name": pl.String,
    }
)

known_reactions_schema = pl.Schema(
        {
        "id": pl.String,
        "smarts": pl.String,
        "enzymes": pl.List(pl.String),
        "reverse": pl.String,
        "db_ids": pl.List(pl.String),
    }
)