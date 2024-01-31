%load_ext autoreload
%autoreload 2

from src.rxn_ctr_mcs import *
from src.utils import load_json, rxn_entry_to_smarts, rm_atom_map_num
from src.pathway_utils import get_reverse_paths_to_starting, create_graph_from_pickaxe
from src.post_processing import *

from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash

from rdkit.Chem import AllChem

from collections import defaultdict
import pandas as pd
import csv
import pickle
import subprocess

# Load processed expansion
starters = 'ccm_v0'
targets = 'mvacid'
generations = 4

expansion_dir = '../data/processed_expansions/'
thermo_dir = '../data/thermo/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name
rxns_path = expansion_dir + 'predicted_reactions_' + fn
paths_path = expansion_dir + 'paths_' + fn
pruned_dir = '../data/pruned_expansions/'

# Load reactions and paths
with open(rxns_path, 'rb') as f:
    pred_rxns = pickle.load(f)

with open(paths_path, 'rb') as f:
    paths = pickle.load(f)

# Load pruned expansion object
pk = Pickaxe()
path = pruned_dir + fn
pk.load_pickled_pickaxe(path)

for k,v in paths.items():
    print(k, len(v))