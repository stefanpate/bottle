import os
import streamlit as st
import polars as pl
from pathlib import Path
from src.post_processing import PathWrangler
from path_viewer.components import HASH_UB, get_existing_usernames, load_user_feedback

CASP_STUDY_ROOT = Path(os.environ.get("CASP_STUDY_ROOT", "/data/processed"))
KNOWN_BIOCHEM_ROOT = Path(os.environ.get("KNOWN_BIOCHEM_ROOT", "/data/known"))

if not CASP_STUDY_ROOT.is_dir():
    st.error(f"CASP_STUDY_ROOT does not exist: {CASP_STUDY_ROOT}")
    st.stop()
if not KNOWN_BIOCHEM_ROOT.is_dir():
    st.error(f"KNOWN_BIOCHEM_ROOT does not exist: {KNOWN_BIOCHEM_ROOT}")
    st.stop()

st.set_page_config(layout="wide")

available_studies = sorted(
    p.name for p in CASP_STUDY_ROOT.iterdir() if p.is_dir()
)

st.sidebar.markdown("### Select CASP Study")
if not available_studies:
    st.sidebar.error("No studies found")
    st.stop()

selected_study = st.sidebar.selectbox(
    "Study",
    options=available_studies,
    index=(
        available_studies.index(st.session_state["study"].name)
        if "study" in st.session_state and st.session_state["study"].name in available_studies
        else 0
    ),
    key="_study_select",
    label_visibility="collapsed",
)

st.sidebar.markdown("---")

study = CASP_STUDY_ROOT / selected_study
known = KNOWN_BIOCHEM_ROOT

study_changed = (
    "study" not in st.session_state
    or st.session_state["study"] != study
)
if study_changed:
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

    st.session_state["study"] = study
    st.session_state["known"] = known
    st.session_state["pw"] = PathWrangler(study=study, known=known)
    st.session_state["username"] = "guest"
    st.session_state["path_feedback"] = {}
    st.session_state["pred_rxn_feedback"] = {}

# Navigation
pages = [
    st.Page("pages/list_view.py", title="Path Viewer", default=True),
    st.Page("pages/path_feedback_review.py", title="Path Feedback Review"),
    st.Page("pages/rxn_feedback_review.py", title="Reaction Feedback Review"),
]

nav = st.navigation(pages)
nav.run()

# Login UI (sidebar, below Apply button added by page scripts)
existing_users = get_existing_usernames(study)

st.sidebar.markdown("---")
st.sidebar.markdown("### Login")

if existing_users:
    selected_user = st.sidebar.selectbox(
        "Select existing user",
        options=existing_users,
        index=None,
        placeholder="Choose a user...",
        key="_login_select",
    )
    if st.sidebar.button("Log in", key="login_select_btn"):
        if selected_user:
            st.session_state["username"] = selected_user
            path_fb, rxn_fb = load_user_feedback(selected_user, study)
            st.session_state["path_feedback"] = path_fb
            st.session_state["pred_rxn_feedback"] = rxn_fb
            st.rerun()
else:
    st.sidebar.info("No existing users found.")

new_user = st.sidebar.text_input("Or enter new username", key="_login_new_user")
if st.sidebar.button("Log in", key="login_new_btn"):
    if new_user and new_user.strip():
        username = new_user.strip()
        st.session_state["username"] = username
        path_fb, rxn_fb = load_user_feedback(username, study)
        st.session_state["path_feedback"] = path_fb
        st.session_state["pred_rxn_feedback"] = rxn_fb
        st.rerun()

current_user = st.session_state["username"]
if current_user == "guest":
    st.sidebar.warning("Browsing as guest (feedback is session-only)")
else:
    st.sidebar.success(f"Logged in as **{current_user}**")
    if st.sidebar.button("Log out"):
        st.session_state["username"] = "guest"
        st.session_state["path_feedback"] = {}
        st.session_state["pred_rxn_feedback"] = {}
        st.rerun()

liked_paths = [pid for pid, v in st.session_state['path_feedback'].items() if v == 1]
if liked_paths:
    st.sidebar.markdown("### My Paths")
    st.sidebar.dataframe(
        pl.DataFrame({"Path ID": [pid[:HASH_UB] for pid in liked_paths]}),
        hide_index=True,
        use_container_width=True,
    )
