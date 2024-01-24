from src.pathway_utils import get_reverse_paths_to_starting, create_graph_from_pickaxe
from src.post_processing import *
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash
from collections import defaultdict

# Set params

expansion_dir = '../data/raw_expansions/'
starters = 'ccm_v0'
targets = 'methylene_molecules'
generations = 4
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name

# Load raw expansion object
pk = Pickaxe()
path = expansion_dir + fn
pk.load_pickled_pickaxe(path)

# # Create the initial graph

# DG, rxn, edge = create_graph_from_pickaxe(pk, "Biology")
# starting_nodes = []
# bad_nodes = []
# for n in DG.nodes():
#     try:
#         if DG.nodes()[n]["Type"] == "Starting Compound":
#             starting_nodes.append(n)
#     except:
#         bad_nodes.append(n)

# # Get pathways
# max_depth = generations * 2
# paths = defaultdict(list)

# Specify Targets / Starting Cpds
# target_cids = [get_compound_hash(smi)[0] for smi in pk.target_smiles]
# target_names = [target_smi_name.loc[smi, "id"] for smi in pk.target_smiles]
target_cids, target_names = [], []
for k,v in pk.targets.items():
    target_cids.append(get_compound_hash(v['SMILES'])[0])
    target_names.append(v['ID'])

starting_cpds = [get_compound_hash(val["SMILES"])[0] for val in pk.compounds.values() if val["Type"].startswith("Start")]

# Loop through targets and get pathways from targets to starting compounds
for i, this_target in enumerate(target_cids):
    this_paths = get_reverse_paths_to_starting(DG, begin_node=this_target, end_nodes=starting_cpds, max_depth=max_depth)
    # If we find paths then reverse those paths and assign to a dictionary
    if this_paths:
        this_paths = list(set([tuple(path[1::2]) for path in [[*reversed(ind_path)] for ind_path in this_paths]]))
        for elt in this_paths:
            for r in pk.reactions[elt[0]]["Reactants"]:
                if r[-1] in starting_cpds:
                    s_name = pk.compounds[r[-1]]["ID"]
                    t_name = target_names[i]
                    paths[(s_name, t_name)].append(pathway(rhashes=elt, starter_hash=r[-1], target_hash=this_target)) 