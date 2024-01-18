from multiprocessing import Process, Queue, Pool
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

# Load reactions and paths
with open(rxns_path, 'rb') as f:
    pred_rxns = pickle.load(f)

# initialize multiprocessing Queue
qin = Queue()
qout = Queue()

for k, v in pred_rxns.items():
    qin.put((k, v))

workers = [Process(target=thermo.rxn_thermo_calcs, args=(qin, qout)) for i in range(n_workers)]

for elt in workers:
    elt.start()

# for elt in workers:
#     elt.join()

ctr = 0
while not qout.empty():
    rxn_id, rxn = qout.get()
    pred_rxns[rxn_id] = rxn
    ctr += 1

assert ctr == len(pred_rxns)

# with open(rxns_path, 'wb') as f:
#     pickle.dump(pred_rxns, f)
    
with open('../data/test_dgrs.pk', 'wb') as f:
    pickle.dump(pred_rxns, f)

print('saved')




