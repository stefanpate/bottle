# import requests

# with open("../artifacts/rheas_to_get_names.txt", 'r') as f:
#     rheas_to_get = [line.strip() for line in f]

# def get_from_mol_elt(mol_elt, keys=('smiles', 'label')):
#     return {k: mol_elt[k] for k in keys}

# def get_from_side(side, keys=('smiles', 'label')):
#     mol_pulls = []
#     for elt in side:
#         if 'reactivePart' in elt:
#             for mol_elt in elt['reactivePart']:
#                 mol_pulls.append(get_from_mol_elt(mol_elt, keys))
#         else:
#             mol_pulls.append(get_from_mol_elt(elt, keys))
#     return mol_pulls

# def get_from_rxn_response(response, keys=('smiles', 'label')):
#     pull = {}
#     pull['left'] = get_from_side(response['left'], keys)
#     pull['right'] = get_from_side(response['right'], keys)
#     return pull


# for i, rid in enumerate(rheas_to_get):
#     url = f"https://www.rhea-db.org/rhea/{rid}/json"
#     response = requests.get(url).json()
#     print(get_from_rxn_response(response))


import requests
import multiprocessing as mp
import json
import time

with open("../artifacts/rheas_to_get_names.txt", 'r') as f:
    rheas_to_get = [line.strip() for line in f]

def get_from_mol_elt(mol_elt, keys=('smiles', 'label')):
    return {k: mol_elt[k] for k in keys}

def get_from_side(side, keys=('smiles', 'label')):
    mol_pulls = []
    for elt in side:
        if 'reactivePart' in elt:
            for mol_elt in elt['reactivePart']:
                mol_pulls.append(get_from_mol_elt(mol_elt, keys))
        else:
            mol_pulls.append(get_from_mol_elt(elt, keys))
    return mol_pulls

def get_from_rxn_rid(rid, keys=('smiles', 'label')):
    try:
        url = f"https://www.rhea-db.org/rhea/{rid}/json"
        response = requests.get(url).json()
        pull = {}
        pull['left'] = get_from_side(response['left'], keys)
        pull['right'] = get_from_side(response['right'], keys)
        return pull
    except:
        return None

n = 100
rheas_to_get = rheas_to_get[:n]

tic = time.perf_counter()
with mp.Pool() as pool:
    res = pool.map(get_from_rxn_rid, rheas_to_get)

toc = time.perf_counter()
print(f"{n} in {toc - tic}s")

res = dict(zip(rheas_to_get, res))

with open("rhea_smiles_names.json", 'w') as f:
    json.dump(res, f)