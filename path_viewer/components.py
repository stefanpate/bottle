import streamlit as st
import polars as pl
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from src.chem_draw import draw_reaction
from src.schemas import path_feedback_schema, rxn_feedback_schema

HASH_UB = 7


def get_existing_usernames(study: Path) -> list[str]:
    usernames = set()
    for fname in ["path_feedback.parquet", "reaction_feedback.parquet"]:
        p = study / fname
        if p.exists():
            df = pl.read_parquet(p, columns=["username"])
            if "username" in df.columns:
                usernames.update(df["username"].unique().to_list())
    return sorted(usernames)


def load_user_feedback(username: str, study: Path) -> tuple[dict, dict]:
    path_fb, rxn_fb = {}, {}
    if username == "guest":
        return path_fb, rxn_fb
    for fname, fb_dict in [("path_feedback.parquet", path_fb), ("reaction_feedback.parquet", rxn_fb)]:
        p = study / fname
        if p.exists():
            df = pl.read_parquet(p)
            if "username" in df.columns:
                df = df.filter(pl.col("username") == username)
                fb_dict.update(dict(zip(df["id"].to_list(), df["feedback"].to_list())))
    return path_fb, rxn_fb


def save_feedback(feedback_dict: dict, filepath: Path, schema: pl.Schema, username: str):
    if username == "guest":
        return
    now = datetime.now()
    new_rows = pl.DataFrame(
        [{"username": username, "id": k, "feedback": v, "date": now.date(), "time": now.time()} for k, v in feedback_dict.items()],
        schema=schema,
    )
    if filepath.exists():
        existing = pl.read_parquet(filepath)
        if "username" not in existing.columns:
            merged = new_rows
        else:
            merged = pl.concat([existing, new_rows]).unique(subset=["username", "id"], keep="last")
    else:
        merged = new_rows
    merged.write_parquet(filepath)


def store_value(key):
    st.session_state[key] = st.session_state["_" + key]


def display_path_metrics(path_id: str, paths_df: pl.DataFrame):
    path_metrics = paths_df.filter(
        pl.col("id") == path_id
    ).select(
        pl.col("mdf"),
        pl.col("mean_max_rxn_sim"),
        pl.col("min_max_rxn_sim"),
        pl.col("feasibility_frac")
    ).to_dicts()[0]

    col1, col2, col3, col4 = st.columns(4)
    mdf = f"{path_metrics['mdf']:.2f}" if path_metrics['mdf'] is not None else "N/A"
    mean_max_sim = f"{path_metrics['mean_max_rxn_sim']:.2f}" if path_metrics['mean_max_rxn_sim'] is not None else "N/A"
    min_max_sim = f"{path_metrics['min_max_rxn_sim']:.2f}" if path_metrics['min_max_rxn_sim'] is not None else "N/A"
    feasibility_frac = f"{path_metrics['feasibility_frac']:.2f}" if path_metrics['feasibility_frac'] is not None else "N/A"
    with col1:
        st.metric("Max-min Driving Force (kJ/mol)", mdf)
    with col2:
        st.metric("Mean Reaction Similarity", mean_max_sim)
    with col3:
        st.metric("Min Reaction Similarity", min_max_sim)
    with col4:
        st.metric("Feasibility Fraction", feasibility_frac)


def display_overall_reaction(prids: list[str], predicted_reactions_smarts: dict[str, str]):
    rxns = [predicted_reactions_smarts[prid] for prid in prids if prid in predicted_reactions_smarts]
    overall_stoich = defaultdict(int)
    for rxn in rxns:
        lhs, rhs = [side.split(".") for side in rxn.split(">>")]
        for rct in lhs:
            overall_stoich[rct] -= 1
        for prd in rhs:
            overall_stoich[prd] += 1

    overall_lhs = []
    overall_rhs = []
    for smi, stoich in overall_stoich.items():
        if stoich < 0:
            for _ in range(-stoich):
                overall_lhs.append(smi)
        elif stoich > 0:
            for _ in range(stoich):
                overall_rhs.append(smi)

    overall_rxn = ".".join(overall_lhs) + ">>" + ".".join(overall_rhs)
    orxn = draw_reaction(overall_rxn).to_str().decode("utf-8")
    st.image(orxn)


def display_predicted_reaction(i: int, prid: str, smarts: str):
    st.write(f"Predicted Reaction {i+1} ({prid[:HASH_UB]})")
    svg = draw_reaction(smarts).to_str().decode("utf-8")
    st.image(svg)


def display_analogue(i: int, ks: pl.DataFrame):
    analogue_select = st.selectbox(
        "Select Analogue",
        key=f"_analogue_select_{i}",
        options=ks["analogue_ids"].to_list(),
        index=0,
        on_change=store_value,
        args=(f"analogue_select_{i}",),
        format_func=lambda x: x[:HASH_UB] if x else "No analogues"
    )
    if analogue_select is not None:
        row = ks.filter(pl.col("analogue_ids") == analogue_select)
        st.write(f"Similarity: {row['rxn_sims'].item():.3f}")
        svg = draw_reaction(row["smarts"].item()).to_str().decode("utf-8")
        st.image(svg)


def display_enzymes(i: int, enz: pl.DataFrame):
    st.dataframe(enz)


def get_path_snapshot(path_id: str, paths_df: pl.DataFrame, predicted_reactions_df: pl.DataFrame, known_reactions_df: pl.DataFrame, enzymes_df: pl.DataFrame):
    if not path_id:
        return None

    prids = paths_df.filter(
        pl.col("id") == path_id
    ).sort(
        by="generation",
        descending=False
    )["rxn_id"].to_list()

    krids_sims = [
        predicted_reactions_df.filter(
            pl.col("id") == prid
        ).select(
            pl.col("analogue_ids").explode(),
            pl.col("rxn_sims").explode()
        ).join(
            known_reactions_df.select(pl.col("id"), pl.col("smarts")),
            left_on="analogue_ids",
            right_on="id",
            how="left",
        ) for prid in prids
    ]

    enz_ids = [
        known_reactions_df.filter(
            pl.col("id").is_in(ks["analogue_ids"].to_list())
        )["enzymes"].explode().unique().to_list()
        for ks in krids_sims
    ]

    enzymes = [
        enzymes_df.filter(
            pl.col("id").is_in(ei)
        ) for ei in enz_ids
    ]
    return prids, krids_sims, enzymes
