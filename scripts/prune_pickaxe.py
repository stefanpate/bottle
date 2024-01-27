from src.pathway_utils import get_reverse_paths_to_starting, create_graph_from_pickaxe
from src.post_processing import *
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash
from collections import defaultdict
import pickle

# Set params
starters = 'methylene_molecules'
targets = 'mvacid'
generations = 4
expansion_dir = '../data/raw_expansions/'
pruned_dir = '../data/pruned_expansions/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name

# Load raw expansion object
pk = Pickaxe()
path = expansion_dir + fn
pk.load_pickled_pickaxe(path)

# Create the initial graph

DG, rxn, edge = create_graph_from_pickaxe(pk, "Biology")
starting_nodes = []
bad_nodes = []
for n in DG.nodes():
    try:
        if DG.nodes()[n]["Type"] == "Starting Compound":
            starting_nodes.append(n)
    except:
        bad_nodes.append(n)

# Get pathways
max_depth = generations * 2
paths = defaultdict(list)

# Specify Targets / Starting Cpds
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

# For pruned pk obj
pruned_rxns = set()
pruned_cpds = set()

# Make predicted reaction dict for processed expansion
pred_rxns = {}
for st_pair in paths:
    for elt in paths[st_pair]:
        for this_rhash in elt.rhashes:
            pruned_rxns.add(this_rhash)

            if this_rhash not in pred_rxns:
                rxn_sma = rxn_hash_2_rxn_sma(this_rhash, pk)
                smi2pkid = get_smi2pkid(this_rhash, pk)
                pred_rxns[this_rhash] = reaction(this_rhash, rxn_sma, smi2pkid)

                for v in smi2pkid.values():
                    pruned_cpds.add(v)

# Save pruned pks
pk.reactions = {k:pk.reactions[k] for k in pruned_rxns}
pk.compounds = {k:pk.compounds[k] for k in pruned_cpds}
pk.pickle_pickaxe(pruned_dir + fn)


# Save reactions dict and paths list (ultimately will replace with expansion object)
rxns_fn = 'predicted_reactions_' + fn
paths_fn = 'paths_' + fn
save_dir = '../data/processed_expansions/'
rxns_path = save_dir + rxns_fn
paths_path = save_dir + paths_fn

with open(rxns_path, 'wb') as f:
    pickle.dump(pred_rxns, f)

with open(paths_path, 'wb') as f:
    pickle.dump(paths, f)

print(f"Saved {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
