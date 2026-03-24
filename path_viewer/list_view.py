import streamlit as st
from pathlib import Path
from src.post_processing import PathWrangler
import polars as pl
import argparse
from src.chem_draw import draw_reaction
from collections import defaultdict

parser = argparse.ArgumentParser(description="Process CASP study directory name.")
parser.add_argument("--casp-study", required=True, help="Name of the CASP study directory")
args = parser.parse_args()

st.set_page_config(layout="wide")
HASH_UB = 7  # Upper bound for displaying hash prefixes

study = Path(f"/home/stef/quest_data/bottle/data/processed/{args.casp_study}")
known = Path("/home/stef/bottle/artifacts/known")

# Callbacks

def get_path_tables():
    path_tables = st.session_state.pw.get_paths(
        starters=st.session_state.starters,
        targets=st.session_state.targets,
        filter_by_enzymes={'existence': st.session_state['loe']},
        sort_by=st.session_state['sort_by'],
    )
    for k, v in path_tables.items():
        st.session_state[k] = v
    
    st.session_state['path_ids'] = st.session_state.paths.unique(subset=['id'], maintain_order=True)['id'].to_list()

@st.cache_data
def get_path_snapshot(path_id: str) -> tuple[list[str], list[pl.DataFrame], list[pl.DataFrame]] | None:
    if not path_id:
        return None
    else:
        prids = st.session_state.paths.filter(
            pl.col("id") == path_id
        ).sort(
            by="generation",
            descending=False
        )["rxn_id"].to_list()

        krids_sims = [
            st.session_state.predicted_reactions.filter(
                pl.col("id") == prid
            ).select(
                pl.col("analogue_ids").explode(),
                pl.col("rxn_sims").explode()
            ) for prid in prids
        ]

        enz_ids = [
            st.session_state.known_reactions.filter(
                pl.col("id").is_in(ks["analogue_ids"].to_list())
            )["enzymes"].explode().unique().to_list() 
            for ks in krids_sims
        ]

        enzymes = [
            st.session_state.enzymes.filter(
            pl.col("id").is_in(ei)
            ) for ei in enz_ids
        ]
        return prids, krids_sims, enzymes


@st.cache_data
def display_predicted_reaction(i: int, prid: str):
    st.write(f"Predicted Reaction {i+1}: {prid[:HASH_UB]}")
    st.image(study / "svgs" / f"{prid}.svg")

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
        sim_value = ks.filter(pl.col("analogue_ids") == analogue_select)["rxn_sims"]
        st.write(f"Similarity: {sim_value.item():.3f}")
        st.image(study / "svgs" / f"{analogue_select}.svg")

def display_enzymes(i: int, enz: pl.DataFrame):
    st.dataframe(enz)

def display_path_metrics(path_id: str):
    path_metrics = st.session_state.paths.filter(
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

@st.cache_data
def display_overall_reaction(prids: list[str]):
    rxns = st.session_state.predicted_reactions.filter(
        pl.col("id").is_in(prids)
    )['smarts'].to_list()
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
            for i in range(-stoich):
                overall_lhs.append(smi)
        elif stoich > 0:
            for i in range(stoich):
                overall_rhs.append(smi)

    overall_rxn = ".".join(overall_lhs) + ">>" + ".".join(overall_rhs)

    fake = draw_reaction(overall_rxn).to_str().decode("utf-8")
    st.image(fake)

def handle_user_input():
    get_path_tables()

def store_value(key):
    st.session_state[key] = st.session_state["_" + key]

def load_value(key):
    st.session_state["_" + key] = st.session_state[key]

def reinitialize_session_state(defaults: dict[str, any]):
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v
        
        load_value(k)

# Initialize session state

pw = PathWrangler(study=study, known=known)

sort_by_options = {
    "mean_max_rxn_sim": "Mean Reaction Similarity",
    "min_max_rxn_sim": "Min Reaction Similarity",
    "mdf": "Max-min Driving Force",
    "feasibility_frac": "Feasibility Fraction"
}

enz_loe_options = [elt.value for elt in pw.enzyme_existence]

defaults = {
    "pw": pw,
    "starters": pw.starters,
    "targets": pw.targets,
    "sort_by": "mean_max_rxn_sim",
    "loe": enz_loe_options[:2],
    "selected_path": None,
}

reinitialize_session_state(defaults)

# UI

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

# Data display panel

selected_path = st.selectbox(
    "Select Path",
    options=st.session_state.get('path_ids', []),
    index=0 if st.session_state.get('path_ids', []) else None,
    key="_selected_path",
    on_change=store_value,
    args=("selected_path",),
    format_func=lambda x: x[:HASH_UB]
)

if st.session_state.selected_path:
    prids, krids_sims, enzymes = get_path_snapshot(st.session_state.selected_path)
    display_path_metrics(st.session_state.selected_path)
    st.header("Overall Reaction")
    display_overall_reaction(prids)

if st.session_state.selected_path:
    header_left, header_right = st.columns([0.6, 0.4])
    with header_left:
        st.header("Predicted Reactions")
        st.write(st.session_state.selected_path)
    with header_right:
        st.header("Known Analogues")

    for i, (prid, ks, enz) in enumerate(zip(prids, krids_sims, enzymes)):
        col_left, col_right = st.columns([0.6, 0.4], border=True)
        with col_left:
            display_predicted_reaction(i, prid)
        with col_right:
            tab_analogue, tab_enzyme = st.tabs(["Analogues", "Enzymes"])
            with tab_analogue:
                display_analogue(i, ks)
            with tab_enzyme:
                display_enzymes(i, enz)
else:
    st.write("No paths found with current criteria.")