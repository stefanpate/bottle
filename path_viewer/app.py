import logging
import os
import streamlit as st
import polars as pl
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s %(message)s",
)
from src.post_processing import PathWrangler
from path_viewer.components import HASH_UB, load_user_feedback, render_auth_header, store_value
from path_viewer.backup_scheduler import start_backup_thread_once

CASP_STUDY_ROOT = Path(os.environ.get("CASP_STUDY_ROOT", "/data/processed"))
KNOWN_BIOCHEM_ROOT = Path(os.environ.get("KNOWN_BIOCHEM_ROOT", "/data/known"))

if not CASP_STUDY_ROOT.is_dir():
    st.error(f"CASP_STUDY_ROOT does not exist: {CASP_STUDY_ROOT}")
    st.stop()
if not KNOWN_BIOCHEM_ROOT.is_dir():
    st.error(f"KNOWN_BIOCHEM_ROOT does not exist: {KNOWN_BIOCHEM_ROOT}")
    st.stop()

st.set_page_config(layout="wide")

start_backup_thread_once()

available_studies = sorted(
    p.name for p in CASP_STUDY_ROOT.iterdir() if p.is_dir()
)

st.sidebar.markdown("### Select CASP Study")
if not available_studies:
    st.sidebar.error("No studies found")
    st.stop()

if "study_select" not in st.session_state:
    st.session_state["study_select"] = available_studies[0] if available_studies else None

if "_study_select" not in st.session_state:
    st.session_state["_study_select"] = available_studies[0] if available_studies else None

st.sidebar.selectbox(
    "Study",
    options=available_studies,
    key="_study_select",
    label_visibility="collapsed",
    on_change=store_value,
    args=("study_select",),
)

st.sidebar.markdown("---")

selected_study = st.session_state["study_select"]
study = CASP_STUDY_ROOT / selected_study
known = KNOWN_BIOCHEM_ROOT

st.session_state["casp_study_root"] = CASP_STUDY_ROOT
st.session_state["known"] = known

username = st.user.email if st.user.is_logged_in else "guest"
st.session_state["username"] = username

study_changed = st.session_state.get("_pw_study") != selected_study
auth_changed = st.session_state.get("_auth_user") != username
if study_changed or auth_changed:
    for k in (
        "starters", "_starters",
        "targets", "_targets",
        "sort_by", "_sort_by",
        "loe", "_loe",
        "mdf_lb", "_mdf_lb",
        "paths", "predicted_reactions", "known_reactions", "enzymes",
        "path_ids", "selected_path", "_selected_path",
    ):
        st.session_state.pop(k, None)

    st.session_state["_pw_study"] = selected_study
    st.session_state["_auth_user"] = username
    st.session_state["pw"] = PathWrangler(study=study, known=known)
    path_fb, rxn_fb = load_user_feedback(username, selected_study)
    st.session_state["path_feedback"] = path_fb
    st.session_state["pred_rxn_feedback"] = rxn_fb

render_auth_header()

# Navigation
pages = [
    st.Page("pages/list_view.py", title="Path Viewer", default=True),
    st.Page("pages/path_feedback_review.py", title="Path Feedback Review"),
    st.Page("pages/rxn_feedback_review.py", title="Reaction Feedback Review"),
]

nav = st.navigation(pages)
nav.run()

liked_paths = [pid for pid, v in st.session_state['path_feedback'].items() if v == 1]
if liked_paths:
    st.sidebar.markdown("### My Paths")
    st.sidebar.dataframe(
        pl.DataFrame({"Path ID": [pid[:HASH_UB] for pid in liked_paths]}),
        hide_index=True,
        use_container_width=True,
    )
