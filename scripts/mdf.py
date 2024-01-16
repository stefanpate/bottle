from multiprocessing import Process, Queue
from src import thermo
import pickle

'''
Params
'''
starters = 'succinate'
targets = 'mvacid'
st_pair = ('succinate', 'mvacid')
generations = 4
num_processes = 1 # Cores to use
################################

expansion_dir = '../data/processed_expansions/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name
rxns_path = expansion_dir + 'predicted_reactions_' + fn
paths_path = expansion_dir + 'paths_' + fn

# Load reactions and paths
with open(rxns_path, 'rb') as f:
    pred_rxns = pickle.load(f)

with open(paths_path, 'rb') as f:
    paths = pickle.load(f)

this_paths = paths[st_pair]

# initialize multiprocessing Queue
q = Queue()

# add all relevant pathway numbers to a shared multiprocessing Queue
for i, elt in enumerate(this_paths):
    q.put(elt)

# define concentration lower and upper bounds (recommended by Dylan Brown of Lucks lab)
lb = 1e-6 # 1 micro-mol
ub = 500e-6 # 500 micro-mol

# start calculating pathway MDF values by multiprocessing
for i in range(num_processes):

    p = Process(target = thermo.calc_pathway_mdf, args=(q,
                                                        lb,
                                                        ub,
                                                        pred_rxns,
                                                        )).start()




