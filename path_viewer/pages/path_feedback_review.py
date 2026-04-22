import streamlit as st
from path_viewer.components import (
    HASH_UB,
    apply_feedback_change,
    store_value,
    display_path_metrics,
    display_overall_reaction,
    display_predicted_reaction,
    display_analogue,
    display_enzymes,
    get_path_snapshot,
)

pw = st.session_state["pw"]
selected_study = st.session_state["study_select"]
study = st.session_state["casp_study_root"] / selected_study
feedback_disabled = not st.user.is_logged_in

st.title("Path Feedback Review")

path_feedback = st.session_state.get("path_feedback", {})

if not path_feedback:
    st.info("No path feedback given yet. Use the Path Viewer to like or dislike paths.")
    st.stop()

# Filter by feedback type
filter_option = st.selectbox(
    "Show",
    options=["All", "Liked", "Disliked"],
    key="path_review_filter",
)

if filter_option == "Liked":
    filtered_ids = [pid for pid, v in path_feedback.items() if v == 1]
elif filter_option == "Disliked":
    filtered_ids = [pid for pid, v in path_feedback.items() if v == 0]
else:
    filtered_ids = list(path_feedback.keys())

if not filtered_ids:
    st.info(f"No {filter_option.lower()} paths found.")
    st.stop()

# Path selector
selected_path = st.selectbox(
    "Select Path",
    options=filtered_ids,
    key="_review_selected_path",
    on_change=store_value,
    args=("review_selected_path",),
    format_func=lambda x: f"{x[:HASH_UB]} ({'liked' if path_feedback.get(x) == 1 else 'disliked'})"
)

if "review_selected_path" not in st.session_state:
    st.session_state["review_selected_path"] = filtered_ids[0] if filtered_ids else None

pid = st.session_state.get("review_selected_path") or selected_path

if pid:
    # Fetch full path data
    path_data = pw.get_path_with_id(pid)
    paths_df = path_data["paths"]
    predicted_reactions_df = path_data["predicted_reactions"]
    known_reactions_df = path_data["known_reactions"]
    enzymes_df = path_data["enzymes"]

    if paths_df.is_empty():
        st.warning(f"Path {pid[:HASH_UB]} not found in data.")
        st.stop()

    snapshot = get_path_snapshot(pid, paths_df, predicted_reactions_df, known_reactions_df, enzymes_df)
    if snapshot:
        prids, krids_sims, enzymes = snapshot

        # Feedback callback for this page
        def store_review_path_feedback(path_id):
            apply_feedback_change(
                st.session_state['path_feedback'],
                path_id,
                st.session_state[f"review_path_fb_{path_id}"],
                study / "path_feedback.parquet",
                st.session_state["username"],
                selected_study,
            )

        def store_review_rxn_feedback(prid):
            apply_feedback_change(
                st.session_state['pred_rxn_feedback'],
                prid,
                st.session_state[f"review_rxn_fb_{prid}"],
                study / "reaction_feedback.parquet",
                st.session_state["username"],
                selected_study,
            )

        # Path metrics
        st.header("Path Metrics")
        display_path_metrics(pid, paths_df)

        # Current feedback + update widget
        if pid in st.session_state['path_feedback']:
            st.session_state[f"review_path_fb_{pid}"] = st.session_state['path_feedback'][pid]

        pred_rxns_smarts = dict(zip(
            predicted_reactions_df["id"].to_list(),
            predicted_reactions_df["smarts"].to_list(),
        ))

        st.header("Overall Reaction")
        display_overall_reaction(prids, pred_rxns_smarts)
        st.caption("Like or dislike this path based on its overall plausibility and usefulness for your goals.")
        st.feedback(
            options="thumbs",
            key=f"review_path_fb_{pid}",
            on_change=store_review_path_feedback,
            args=(pid,),
            disabled=feedback_disabled,
        )
        if feedback_disabled:
            st.caption("Sign in to leave feedback.")

        # Predicted reactions with analogues
        header_left, header_right = st.columns([0.6, 0.4])
        with header_left:
            st.header("Predicted Reactions")
        with header_right:
            st.header("Known Analogues")

        for i, (prid, ks, enz) in enumerate(zip(prids, krids_sims, enzymes)):
            col_left, col_right = st.columns([0.6, 0.4], border=True)
            with col_left:
                display_predicted_reaction(i, prid, pred_rxns_smarts[prid])
                if prid in st.session_state.get('pred_rxn_feedback', {}):
                    st.session_state[f"review_rxn_fb_{prid}"] = st.session_state['pred_rxn_feedback'][prid]
                st.caption("Like or dislike this predicted reaction based on its plausibility.")
                st.feedback(
                    options="thumbs",
                    key=f"review_rxn_fb_{prid}",
                    on_change=store_review_rxn_feedback,
                    args=(prid,),
                    disabled=feedback_disabled,
                )
            with col_right:
                tab_analogue, tab_enzyme = st.tabs(["Analogues", "Enzymes"])
                with tab_analogue:
                    display_analogue(i, ks)
                with tab_enzyme:
                    display_enzymes(i, enz)
