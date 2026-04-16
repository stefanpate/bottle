import streamlit as st
import polars as pl
from src.schemas import rxn_feedback_schema
from path_viewer.components import (
    HASH_UB,
    save_feedback,
    display_predicted_reaction,
    display_analogue,
    display_enzymes,
)

pw = st.session_state["pw"]
study = st.session_state["study"]

st.title("Reaction Feedback Review")

rxn_feedback = st.session_state.get("pred_rxn_feedback", {})

if not rxn_feedback:
    st.info("No reaction feedback given yet. Use the Path Viewer to like or dislike reactions.")
    st.stop()

# Filter by feedback type
filter_option = st.selectbox(
    "Show",
    options=["All", "Liked", "Disliked"],
    key="rxn_review_filter",
)

if filter_option == "Liked":
    filtered_ids = [rid for rid, v in rxn_feedback.items() if v == 1]
elif filter_option == "Disliked":
    filtered_ids = [rid for rid, v in rxn_feedback.items() if v == 0]
else:
    filtered_ids = list(rxn_feedback.keys())

if not filtered_ids:
    st.info(f"No {filter_option.lower()} reactions found.")
    st.stop()

# Get path counts for all filtered reactions
path_counts = pw.count_paths_per_reaction(filtered_ids)

# Reaction selector
selected_rxn = st.selectbox(
    "Select Reaction",
    options=filtered_ids,
    key="review_selected_rxn",
    format_func=lambda x: f"{x[:HASH_UB]} ({'liked' if rxn_feedback.get(x) == 1 else 'disliked'}) - in {path_counts.get(x, 0)} path(s)"
)

if selected_rxn:
    prid = selected_rxn

    # Feedback callback
    def store_review_rxn_feedback(prid):
        st.session_state['pred_rxn_feedback'][prid] = st.session_state[f"rxn_review_fb_{prid}"]
        save_feedback(st.session_state['pred_rxn_feedback'], study / "reaction_feedback.parquet", rxn_feedback_schema, st.session_state["username"])

    # Path count metric
    st.metric("Paths containing this reaction", path_counts.get(prid, 0))

    # Reaction display
    st.header("Predicted Reaction")
    display_predicted_reaction(0, prid, str(study))

    # Feedback widget
    if prid in st.session_state['pred_rxn_feedback']:
        st.session_state[f"rxn_review_fb_{prid}"] = st.session_state['pred_rxn_feedback'][prid]
    st.caption("Like or dislike this predicted reaction based on its plausibility.")
    st.feedback(
        options="thumbs",
        key=f"rxn_review_fb_{prid}",
        on_change=store_review_rxn_feedback,
        args=(prid,),
    )

    # Load analogue and enzyme data for this reaction
    prxn_df = pl.scan_parquet(pw.predicted_reactions).filter(pl.col("id") == prid).collect()

    if not prxn_df.is_empty():
        krids_sims = prxn_df.select(
            pl.col("analogue_ids").explode(),
            pl.col("rxn_sims").explode()
        )

        krxns_df = pl.scan_parquet(pw.known_reactions).filter(
            pl.col("id").is_in(krids_sims["analogue_ids"].to_list())
        ).collect()

        enz_ids = krxns_df["enzymes"].explode().unique().to_list()
        enz_df = pl.scan_parquet(pw.enzymes).filter(pl.col("id").is_in(enz_ids)).collect()

        # Analogues and enzymes
        st.header("Known Analogues")
        tab_analogue, tab_enzyme = st.tabs(["Analogues", "Enzymes"])
        with tab_analogue:
            display_analogue(0, krids_sims, str(study))
        with tab_enzyme:
            display_enzymes(0, enz_df)
