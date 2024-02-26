'''
1. Find paths between starters and targets
2. Create 

'''
from src.pathway_utils import get_reverse_paths_to_starting, create_graph_from_pickaxe
from src.utils import load_json
from src.post_processing import *
from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash
from collections import defaultdict
import pickle
import csv
import pandas as pd

# Set params
starters = '2mg'
targets = 'mvacid'
generations = 2
expansion_dir = '../data/raw_expansions/' # To load
pruned_dir = '../data/pruned_expansions/' # To save
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000" # Expansion file name

# Load raw expansion object
pk = Pickaxe()
path = expansion_dir + fn + '.pk'
pk.load_pickled_pickaxe(path)

pe = ProcessedExpansion() # Initialize processed expansion

# For pruned pk obj
pruned_rxns = set()
pruned_cpds = set()

'''
Load known reaction information
'''
known_rxns = load_json("../data/mapping/known_rxns_w_provenance_all_info_jni.json")
rule2kr = defaultdict(lambda : defaultdict(lambda : {'db_ids':[], 'smarts':None, 'uniprot_ids':[], 'imt_rules':[]})) # Look up known rxns by rule
for kr, db_entries in known_rxns.items():
    for db_id, entry in db_entries.items():
        for rule in entry['imt_rules']:
            rule2kr[rule][kr]['db_ids'].append(db_id)
            rule2kr[rule][kr]['uniprot_ids'] += entry['uniprot_ids']
            rule2kr[rule][kr]['imt_rules'] = entry['imt_rules'] # Confusing

            if rule2kr[rule][kr]['smarts'] is None:
                rule2kr[rule][kr]['smarts'] = entry['smarts']


'''
Get pathways to target molecule
'''


def construct_pr_list(pk_rids, rule2kr, pk):
    prs = []
    for rid in pk_rids:
        pk_reaction = pk.reactions[rid]

        # Collect known reaction analogues inffo
        krs = []
        for rule in pk_reaction["Operators"]:
            for known_rid, entry in rule2kr[rule].items():
                db_entries = [DatabaseEntry(x.split(':')[0], x.split(":")[1]) for x in entry['db_ids']]
                enzymes = [Enzyme(x, None) for x in entry["uniprot_ids"]]
                kr = KnownReaction(id=known_rid, smarts=entry['smarts'],
                    imt_rules=entry['imt_rules'], database_entries=db_entries,
                    enzymes=enzymes)
                krs.append(kr)
        
        # Construct predicted reaction object with analogues
        sma = rxn_hash_2_rxn_sma(rid, pk)
        pr = PredictedReaction(id=rid, smarts=sma, imt_rules=list(pk_reaction["Operators"]), analogues=krs)
        prs.append(pr)

    return prs

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
                    prs = construct_pr_list(elt, rule2kr, pk)
                    pe.add_path(s_name, t_name, prs)

                    for pk_rid in elt:
                        pruned_rxns.add(pk_rid)
                        for _, cpd_id in pk.reactions[pk_rid]["Reactants"]:
                            pruned_cpds.add(cpd_id)
                        for _, cpd_id in pk.reactions[pk_rid]["Products"]:
                            pruned_cpds.add(cpd_id)


'''
Save
'''
# Save pruned pks
pk.reactions = {k:pk.reactions[k] for k in pruned_rxns}
pk.compounds = {k:pk.compounds[k] for k in pruned_cpds}
print(f"Pruned pk to {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
print("Saving pruned pk object")
pk.pickle_pickaxe(pruned_dir + fn + '.pk')

# Save processed expansion object
print("Saving processed expansion object")
save_dir = '../data/processed_expansions/'
with open(save_dir + fn + '.pkl', 'wb') as f:
    pickle.dump(pe, f)
