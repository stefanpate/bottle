import os
import sqlite3
import streamlit as st
import polars as pl
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
from src.chem_draw import draw_reaction

HASH_UB = 7

_FILENAME_TO_TABLE = {
    "path_feedback.parquet": "path_feedback",
    "reaction_feedback.parquet": "reaction_feedback",
}


def _db_path() -> Path:
    return Path(os.environ.get("FEEDBACK_ROOT", "/feedback")) / "feedback.db"


def _connect(db_path: Path | None = None) -> sqlite3.Connection:
    path = db_path if db_path is not None else _db_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path, isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA foreign_keys=ON;")
    for table in ("path_feedback", "reaction_feedback"):
        conn.execute(
            f"""CREATE TABLE IF NOT EXISTS {table} (
                username TEXT NOT NULL,
                casp_study TEXT NOT NULL,
                id TEXT NOT NULL,
                feedback INTEGER NOT NULL,
                updated_at TEXT NOT NULL,
                PRIMARY KEY (username, casp_study, id)
            )"""
        )
        conn.execute(
            f"CREATE INDEX IF NOT EXISTS ix_{table}_user_study ON {table}(username, casp_study)"
        )
    _migrate_legacy_parquets(conn)
    return conn


def _migrate_legacy_parquets(conn: sqlite3.Connection) -> None:
    counts = {
        t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        for t in ("path_feedback", "reaction_feedback")
    }
    if any(counts.values()):
        return
    study_root = Path(os.environ.get("CASP_STUDY_ROOT", "/data/processed"))
    if not study_root.is_dir():
        return
    for study_dir in study_root.iterdir():
        if not study_dir.is_dir():
            continue
        for fname, table in _FILENAME_TO_TABLE.items():
            p = study_dir / fname
            if not p.exists():
                continue
            try:
                df = pl.read_parquet(p)
            except Exception:
                continue
            if "username" not in df.columns or "id" not in df.columns or "feedback" not in df.columns:
                continue
            rows = [
                (r["username"], study_dir.name, r["id"], int(r["feedback"]), _iso_from_row(r))
                for r in df.iter_rows(named=True)
            ]
            if rows:
                conn.executemany(
                    f"INSERT OR IGNORE INTO {table} (username, casp_study, id, feedback, updated_at) VALUES (?, ?, ?, ?, ?)",
                    rows,
                )


def _iso_from_row(row: dict) -> str:
    d = row.get("date")
    t = row.get("time")
    if d is not None and t is not None:
        return f"{d}T{t}"
    return datetime.now(timezone.utc).isoformat()


def get_existing_usernames() -> list[str]:
    with _connect() as conn:
        rows = conn.execute(
            "SELECT username FROM path_feedback "
            "UNION SELECT username FROM reaction_feedback"
        ).fetchall()
    return sorted({r[0] for r in rows})


def load_user_feedback(username: str, casp_study: str) -> tuple[dict, dict]:
    if username == "guest":
        return {}, {}
    with _connect() as conn:
        path_rows = conn.execute(
            "SELECT id, feedback FROM path_feedback "
            "WHERE username = ? AND casp_study = ?",
            (username, casp_study),
        ).fetchall()
        rxn_rows = conn.execute(
            "SELECT id, feedback FROM reaction_feedback "
            "WHERE username = ? AND casp_study = ?",
            (username, casp_study),
        ).fetchall()
    return dict(path_rows), dict(rxn_rows)


def save_feedback(feedback_dict: dict, filepath: Path, schema: object, username: str, casp_study: str):
    if username == "guest" or not feedback_dict:
        return
    table = _FILENAME_TO_TABLE.get(Path(filepath).name)
    if table is None:
        raise ValueError(f"Unknown feedback filepath: {filepath}")
    now = datetime.now(timezone.utc).isoformat()
    rows = [(username, casp_study, k, int(v), now) for k, v in feedback_dict.items()]
    with _connect() as conn:
        conn.executemany(
            f"""INSERT INTO {table} (username, casp_study, id, feedback, updated_at)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(username, casp_study, id) DO UPDATE SET
                    feedback = excluded.feedback,
                    updated_at = excluded.updated_at""",
            rows,
        )


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
