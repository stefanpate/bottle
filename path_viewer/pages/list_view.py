import streamlit as st
import polars as pl
from src.schemas import path_feedback_schema, rxn_feedback_schema
from path_viewer.components import (
    HASH_UB,
    save_feedback,
    store_value,
    display_path_metrics,
    display_overall_reaction,
    display_predicted_reaction,
    display_analogue,
    display_enzymes,
    get_path_snapshot,
)

pw = st.session_state["pw"]
study = st.session_state["study"]

sort_by_options = {
    "mean_max_rxn_sim": "Mean Reaction Similarity",
    "min_max_rxn_sim": "Min Reaction Similarity",
    "mdf": "Max-min Driving Force",
    "feasibility_frac": "Feasibility Fraction"
}

enz_loe_options = [elt.value for elt in pw.enzyme_existence]

# Initialize page-specific session state defaults
page_defaults = {
    "starters": pw.starters,
    "targets": pw.targets,
    "sort_by": "mean_max_rxn_sim",
    "loe": enz_loe_options[:2],
}
for k, v in page_defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v
    if "_" + k not in st.session_state:
        st.session_state["_" + k] = st.session_state[k]


# Callbacks

def get_path_tables():
    blacklist_path_ids = [pid for pid, v in st.session_state['path_feedback'].items() if v == 0] or None
    blacklist_rxn_ids = [prid for prid, v in st.session_state['pred_rxn_feedback'].items() if v == 0] or None

    path_tables = pw.get_paths(
        starters=st.session_state.starters,
        targets=st.session_state.targets,
        filter_by_enzymes={'existence': st.session_state['loe']},
        sort_by=st.session_state['sort_by'],
        blacklist_path_ids=blacklist_path_ids,
        blacklist_rxn_ids=blacklist_rxn_ids,
    )
    for k, v in path_tables.items():
        st.session_state[k] = v

    st.session_state['path_ids'] = st.session_state.paths.unique(subset=['id'], maintain_order=True)['id'].to_list()
    if st.session_state['path_ids']:
        st.session_state['_selected_path'] = st.session_state['path_ids'][0]
        st.session_state['selected_path'] = st.session_state['path_ids'][0]
    else:
        st.session_state['_selected_path'] = None
        st.session_state['selected_path'] = None


def handle_user_input():
    get_path_tables()


def store_pred_rxn_feedback(prid):
    st.session_state['pred_rxn_feedback'][prid] = st.session_state[f"pred_rxn_feedback_{prid}"]
    save_feedback(st.session_state['pred_rxn_feedback'], study / "reaction_feedback.parquet", rxn_feedback_schema, st.session_state["username"])


def ban_rules(rule_pairs):
    rule_pairs_set = set(rule_pairs)
    prxns = st.session_state.predicted_reactions
    for row in prxns.iter_rows(named=True):
        pairs = set(zip(row["rules"], row["rule_sets"]))
        if pairs & rule_pairs_set:
            st.session_state['pred_rxn_feedback'][row["id"]] = 0
    save_feedback(st.session_state['pred_rxn_feedback'], study / "reaction_feedback.parquet", rxn_feedback_schema, st.session_state["username"])


def store_path_feedback(path_id):
    st.session_state['path_feedback'][path_id] = st.session_state[f"path_feedback_{path_id}"]
    save_feedback(st.session_state['path_feedback'], study / "path_feedback.parquet", path_feedback_schema, st.session_state["username"])


# Sidebar

starters = st.sidebar.multiselect(
    "Starting Compounds",
    options=st.session_state.starters,
    key="_starters",
    on_change=store_value,
    args=("starters",)
)

targets = st.sidebar.multiselect(
    "Target Compounds",
    options=st.session_state.targets,
    key="_targets",
    on_change=store_value,
    args=("targets",)
)

sort_by = st.sidebar.selectbox(
    "Sort by",
    options=sort_by_options.keys(),
    format_func=lambda x: sort_by_options[x],
    key="_sort_by",
    on_change=store_value,
    args=("sort_by",)
)

loe = st.sidebar.multiselect(
    "Enzyme Level of Evidence",
    options=enz_loe_options,
    key="_loe",
    on_change=store_value,
    args=("loe",)
)

apply = st.sidebar.button(
    "Apply",
    on_click=handle_user_input
)

liked_paths = [pid for pid, v in st.session_state['path_feedback'].items() if v == 1]
if liked_paths:
    st.sidebar.markdown("### My Paths")
    st.sidebar.dataframe(
        pl.DataFrame({"Path ID": [pid[:HASH_UB] for pid in liked_paths]}),
        hide_index=True,
        use_container_width=True,
    )

# Data display panel
selected_path = st.selectbox(
    "Select Path",
    options=st.session_state.get('path_ids', []),
    key="_selected_path",
    on_change=store_value,
    args=("selected_path",),
    format_func=lambda x: x[:HASH_UB]
)

if st.session_state.get('selected_path') and 'paths' in st.session_state:
    snapshot = get_path_snapshot(
        st.session_state.selected_path,
        st.session_state.paths,
        st.session_state.predicted_reactions,
        st.session_state.known_reactions,
        st.session_state.enzymes,
    )
    if snapshot:
        prids, krids_sims, enzymes = snapshot
        st.header("Path Metrics")
        display_path_metrics(st.session_state.selected_path, st.session_state.paths)
        pid = st.session_state.selected_path
        if pid in st.session_state['path_feedback']:
            st.session_state[f"path_feedback_{pid}"] = st.session_state['path_feedback'][pid]

        # Build smarts dict for overall reaction display
        pred_rxns_smarts = dict(zip(
            st.session_state.predicted_reactions["id"].to_list(),
            st.session_state.predicted_reactions["smarts"].to_list(),
        ))

        st.header("Overall Reaction")
        display_overall_reaction(prids, pred_rxns_smarts, str(study))
        st.caption("Like or dislike this path based on its overall plausibility and usefulness for your goals.")
        st.feedback(
            options="thumbs",
            key=f"path_feedback_{pid}",
            on_change=store_path_feedback,
            args=(pid,),
        )

        # Reactions side-by-side with known analogues and enzymes
        header_left, header_right = st.columns([0.6, 0.4])
        with header_left:
            st.header("Predicted Reactions")
        with header_right:
            st.header("Known Analogues")

        for i, (prid, ks, enz) in enumerate(zip(prids, krids_sims, enzymes)):
            col_left, col_right = st.columns([0.6, 0.4], border=True)
            with col_left:
                display_predicted_reaction(i, prid, str(study))
                if prid in st.session_state['pred_rxn_feedback']:
                    st.session_state[f"pred_rxn_feedback_{prid}"] = st.session_state['pred_rxn_feedback'][prid]
                fb_col, ban_col, _ = st.columns([1, 2, 10], vertical_alignment="center")
                st.caption("Like or dislike this predicted reaction based on its plausibility.")
                with fb_col:
                    st.feedback(
                        options="thumbs",
                        key=f"pred_rxn_feedback_{prid}",
                        on_change=store_pred_rxn_feedback,
                        args=(prid,),
                    )
                with ban_col:
                    prxn_row = st.session_state.predicted_reactions.filter(pl.col("id") == prid)
                    rules = prxn_row["rules"].explode().to_list()
                    rule_sets = prxn_row["rule_sets"].explode().to_list()
                    rule_pairs = list(zip(rules, rule_sets))
                    st.button(
                        "Rule",
                        key=f"ban_rule_{prid}",
                        icon=":material/thumb_down:",
                        on_click=ban_rules,
                        args=(rule_pairs,),
                        help="Dislike all reactions generated by this reaction's rule(s)",
                    )
            with col_right:
                tab_analogue, tab_enzyme = st.tabs(["Analogues", "Enzymes"])
                with tab_analogue:
                    display_analogue(i, ks, str(study))
                with tab_enzyme:
                    display_enzymes(i, enz)
elif st.session_state.get('path_ids') is None:
    st.write("Please click 'Apply' to load paths based on selected criteria.")
else:
    st.write("No paths found with current criteria.")
