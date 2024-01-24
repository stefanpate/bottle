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
n_workers = 2 # Cores to use
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
qin = Queue()
qout = Queue()

# add all relevant pathway numbers to a shared multiprocessing Queue
for elt in this_paths[:2]:
    qin.put(elt)


# define concentration lower and upper bounds (recommended by Dylan Brown of Lucks lab)
lb = 1e-6 # 1 micro-mol
ub = 500e-6 # 500 micro-mol

workers = [Process(target=thermo.calc_pathway_mdf, args=(qin, qout, lb, ub, pred_rxns)) for i in range(n_workers)]

for elt in workers:
    elt.start()

for elt in workers:
    elt.join()

new_this_paths = []
while not qout.empty():
    new_this_paths.append(qout.get())

paths[st_pair] = new_this_paths

# with open(rxns_path, 'wb') as f:
#     pickle.dump(pred_rxns, f)

with open('../data/test_thermo_multi.pk', 'wb') as f:
    pickle.dump(paths, f)

print('saved')




