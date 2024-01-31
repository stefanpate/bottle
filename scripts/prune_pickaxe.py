from src.pathway_utils import get_reverse_paths_to_starting, create_graph_from_pickaxe
from src.utils import load_json, rxn_entry_to_smarts
from src.post_processing import *
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash
from collections import defaultdict
import pickle
import csv
import pandas as pd

# Set params
starters = 'ccm_v0'
targets = 'hopa'
generations = 4
expansion_dir = '../data/raw_expansions/'
pruned_dir = '../data/pruned_expansions/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name

# Load raw expansion object
pk = Pickaxe()
path = expansion_dir + fn
pk.load_pickled_pickaxe(path)

'''
Get pathways to target molecule
'''
print("Finding paths")
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

'''
Prune and create processed object
'''
print("Pruning & creating processed object")
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

'''
Map predicted reactions to know reactions
'''
print("Mapping predicted to known reactions")
# Load rules
rules_path = '../data/rules/JN3604IMT_rules.tsv'
rule_df = pd.read_csv(rules_path, delimiter='\t')
rule_df.set_index('Name', inplace=True)

# Load mapping
rxn2rule = {}
db_names = ['_mc_v21', '_brenda', '_kegg']
suffix = '_imt_rules_enforce_cof.csv'
for name in db_names:
    mapping_path = '../data/mapping/mapping' + name + suffix
    with open(mapping_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) == 1:
                rxn2rule[row[0]] = []
            else:
                rxn2rule[row[0]] = row[1:]

# Make rule2rxn
rule2rxn = {}
for k,v in rxn2rule.items():
    for elt in v:
        if elt not in rule2rxn:
            rule2rxn[elt] = [k]
        else:
            rule2rxn[elt].append(k)

# Load all known reaction json entries into dict
known_rxns = {}
pref = '../data/mapping/'
suffs = ['mc_v21_as_is.json', 'brenda_as_is.json', 'kegg_as_is.json']
for elt in suffs:
    known_rxns.update(load_json(pref + elt))

# Populate pred_rxns w/ mapped known reactions
for k, v in pred_rxns.items():
    this_rules = list(pk.reactions[k]["Operators"])
    this_known_rxns = []
    for elt in this_rules:
        if elt in rule2rxn:
            this_rxn_ids = rule2rxn[elt]
            for this_id in this_rxn_ids:
                this_sma = rxn_entry_to_smarts(known_rxns[this_id])
                this_known_rxns.append((None, this_sma, this_id))
    
    v.known_rxns = [list(elt) for elt in set(this_known_rxns)]

'''
Save
'''
# Save pruned pks
print("Saving pruned pk object")
pk.reactions = {k:pk.reactions[k] for k in pruned_rxns}
pk.compounds = {k:pk.compounds[k] for k in pruned_cpds}
pk.pickle_pickaxe(pruned_dir + fn)

# Save reactions dict and paths list (ultimately will replace with expansion object)
print("Saving processed expansion object")
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
