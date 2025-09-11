import streamlit as st
from pathlib import Path
import os
from src.post_processing import PathWrangler

# Initialize path wrangler
study = Path("/home/stef/quest_data/bottle/data/processed/test")
known = Path("/home/stef/bottle/artifacts/known")
pw = PathWrangler(study=study, known=known)
st.session_state['pw'] = pw

# Callbacks

@st.cache_data
def get_path_tables():
    starter_select = [st.session_state.starter_lookup[elt] for elt in st.session_state.starters]
    target_select = [st.session_state.target_lookup[elt] for elt in st.session_state.targets]
    path_tables = st.session_state.pw.get_paths(
        starters=starter_select,
        targets=target_select,
        filter_by_enzymes={'existence': st.session_state['loe']},
        sort_by=st.session_state['sort_by'],
    )
    st.session_state['path_tables'] = path_tables

def display_list_view_path():
    pass

# def store_value(key):
#     st.session_state[key] = st.session_state['_' + key]

# def load_value(key):
#     st.session_state['_' + key] = st.session_state[key]

def handle_user_input(key):
    # store_value(key)
    get_path_tables()

# UI

sort_by_options = {
    "mean_max_rxn_sim": "Mean Reaction Similarity",
    "min_max_rxn_sim": "Min Reaction Similarity",
    "mdf": "Max-min Driving Force",
    "feasbility_frac": "Feasibility Fraction"
}

ub = 7
starter_lookup = {f"Starter #{i+1}" : elt for i, elt in enumerate(pw.starters)}
target_lookup = {elt[:ub] : elt for elt in pw.targets}
enz_loe_options = [elt.value for elt in pw.enzyme_existence]

st.session_state['starter_lookup'] = starter_lookup
st.session_state['target_lookup'] = target_lookup

starters = st.sidebar.multiselect(
    "Starting Compounds",
    options=starter_lookup.keys(),
    key="starters",
    on_change=handle_user_input,
    args=("starters",),
    default=starter_lookup.keys()
)

targets = st.sidebar.multiselect(
    "Target Compounds",
    options=target_lookup.keys(),
    key="targets",
    on_change=handle_user_input,
    args=("targets",),
    default=target_lookup.keys()
)

sort_by = st.sidebar.selectbox(
    "Sort by",
    options=sort_by_options.keys(),
    format_func=lambda x: sort_by_options[x],
    key="sort_by"
)

loe = st.sidebar.multiselect(
    "Enzyme Level of Evidence",
    options=enz_loe_options,
    key="loe",
    default=enz_loe_options[:2],
)

# Data display 
col1, col2 = st.columns([0.6, 0.4])

with col1:
    st.header("Predicted Path")
    display_list_view_path()

# print(st.session_state)
# print(st.session_state.get('path_tables', None))
