import streamlit as st
import polars as pl
from pathlib import Path
import argparse
from src.post_processing import PathWrangler
from path_viewer.components import HASH_UB, get_existing_usernames, load_user_feedback

parser = argparse.ArgumentParser(description="Process CASP study directory name.")
parser.add_argument("--casp-study", required=True, help="Name of the CASP study directory")
args = parser.parse_args()

st.set_page_config(layout="wide")

study = Path(f"/home/stef/quest_data/bottle/data/processed/{args.casp_study}")
known = Path("/home/stef/bottle/artifacts/known")

# Initialize core session state
if "pw" not in st.session_state:
    st.session_state["pw"] = PathWrangler(study=study, known=known)
    st.session_state["study"] = study
    st.session_state["known"] = known
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
