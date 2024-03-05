from src.utils import load_json, get_compound_hash, get_reaction_hash, postsanitize_smiles, neutralise_charges
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS
from collections import Counter, defaultdict
from copy import deepcopy
from itertools import permutations
from tqdm import tqdm
import pickle
import logging

def sanitize(reactants, products):
    return postsanitize_smiles(reactants)[0], postsanitize_smiles(products)[0]

def neutralize(reactants, products):
    reactants = [Chem.MolToSmiles(neutralise_charges(Chem.MolFromSmiles(elt, sanitize=False))) for elt in reactants]
    products = [Chem.MolToSmiles(neutralise_charges(Chem.MolFromSmiles(elt, sanitize=False))) for elt in products]
    return reactants, products

def ct_any_ids(known_rxns):
    n_rxns = len(known_rxns)
    n_any_ids = 0
    for k,v in known_rxns.items():
        if v['uniprot_ids']:
            n_any_ids += 1

    return n_any_ids, n_rxns

def smarts_to_sub_smiles(smarts):
    reactants, products = smarts.split(">>")
    reactants = reactants.split('.')
    products = products.split('.')
    return reactants, products

def sub_smiles_to_smarts(reactants, products):
    sma = ".".join(reactants) + '>>' + ".".join(products)
    return sma

def get_rhash(reactants, products):
    reactants_hash_stoich = list(Counter([get_compound_hash(elt)[0] for elt in reactants]).items())
    products_hash_stoich = list(Counter([get_compound_hash(elt)[0] for elt in products]).items())
    rhash = get_reaction_hash(reactants_hash_stoich, products_hash_stoich)
    return rhash

# Load kwown reactions and uniprot entries
known_rxns = load_json("../data/mapping/known_rxns_jni.json")
jni_up = load_json("../data/mapping/jni_uniprot.json")


# Assemble rhea reactions
rhea_smarts = pd.read_csv("../data/mapping/rhea-reaction-smiles.tsv", sep='\t', header=None, names=['rhea_id', 'smarts'])
rhea_directions = pd.read_csv("../data/mapping/rhea-directions.tsv", sep='\t')
rhea_smarts.head()

# See what types of rhea ids I got
rhea_ids_from_up = set()
master_ids = set(rhea_directions.loc[:, "RHEA_ID_MASTER"].values)
lr_ids = set(rhea_directions.loc[:, "RHEA_ID_LR"].values)
rl_ids = set(rhea_directions.loc[:, "RHEA_ID_RL"].values)
bi_ids = set(rhea_directions.loc[:, "RHEA_ID_BI"].values)
ids_w_sma = set(rhea_smarts.loc[:, 'rhea_id'].values)
other_ids = set()

for k,v in known_rxns.items():
    for uid in v['uniprot_ids']:
        if uid in jni_up:
            for rid in jni_up[uid]['rhea_ids']:
                rid = int(rid.split(':')[-1])
                rhea_ids_from_up.add(rid)

                if (rid not in master_ids) and (rid not in lr_ids) and (rid not in rl_ids) and (rid not in bi_ids) and (rid not in ids_w_sma):
                    other_ids.add(rid)

print(f"total needed: {len(rhea_ids_from_up)}, master: {len(master_ids)}, LR: {len(lr_ids)}, RL: {len(rl_ids)}, BI: {len(bi_ids)}, Other: {len(other_ids)}")

'''
For all rows in rhea directions
    Check if lr and rl ids in ids with sma
    If both: get sma, 
    elif one: get sma, reverse it
    else continue
    put in tuple, enter in dict under all ids
'''


# Multiple rhea ids will map to 
# same tuple of reaction info
rhea_id_to_rhashes = {}
rhea_id_to_stoichless_rhashes = {}
rhea_id_to_smarts = {}
for _, row in rhea_directions.iterrows():
    lrid, rlid = row["RHEA_ID_LR"], row["RHEA_ID_RL"]
    if lrid in ids_w_sma and rlid in ids_w_sma:
        lr_sma = rhea_smarts.loc[rhea_smarts['rhea_id'] == lrid, 'smarts'].values[0]
        rl_sma = rhea_smarts.loc[rhea_smarts['rhea_id'] == rlid, 'smarts'].values[0]
    elif lrid in ids_w_sma:
        lr_sma = rhea_smarts.loc[rhea_smarts['rhea_id'] == lrid, 'smarts'].values[0]
        rl_sma = ">>".join([lr_sma.split('>>')[1], lr_sma.split('>>')[0]])
    elif rlid in ids_w_sma:
        rl_sma = rhea_smarts.loc[rhea_smarts['rhea_id'] == rlid, 'smarts'].values[0]
        lr_sma = ">>".join([lr_sma.split('>>')[1], lr_sma.split('>>')[0]])
    else:
        continue

    san_sma = []
    rhashes = []
    for sma in [lr_sma, rl_sma]:
        reactants, products = smarts_to_sub_smiles(sma)
        reactants = [elt for elt in reactants if elt != '[H+]'] # Remove protons
        products = [elt for elt in products if elt != '[H+]']
        reactants, products = sanitize(reactants, products)
        reactants, products = neutralize(reactants, products)
        sma = ".".join(reactants) + '>>' + ".".join(products)
        rhash = get_rhash(reactants, products)
        san_sma.append(sma)
        rhashes.append(rhash)

        for id in row.to_list():
            rhea_id_to_rhashes[id] = tuple(rhashes)
            rhea_id_to_smarts[id] = tuple(san_sma)

'''
Filters
'''
in_uniprot = lambda x : x in jni_up # Whether I found uniprot id in uniprot
has_rxn = lambda x : len(jni_up[x]['rhea_ids']) > 0 # Whether a found uniprot entry had rhea reactions at all

def is_rhash_equivalent(x, putative_rhash, jni_up, rhea_id_to_rhashes):
        '''
        Do any rhea reaction rhashes form up entry match the one from joseph's
        known reactions?
        '''
        rheas = [int(elt.split(":")[1]) for elt in jni_up[x]['rhea_ids']]
        bool_list = [putative_rhash in rhea_id_to_rhashes.get(elt, []) for elt in rheas]
        return any(bool_list)

def is_stoichless_rhash_equivalent(x, putative_rhash, jni_up, rhea_id_to_rhashes):
        '''
        Do any rhea reaction rhashes form up entry match the one from joseph's
        known reactions?
        '''
        rheas = [int(elt.split(":")[1]) for elt in jni_up[x]['rhea_ids']]
        bool_list = [putative_rhash in rhea_id_to_rhashes.get(elt, []) for elt in rheas]
        return any(bool_list)

def mcs_index(smi1, smi2, do_valence, norm):
    mol1, mol2 = [Chem.MolFromSmiles(elt, sanitize=True) for elt in [smi1, smi2]]
    
    if mol1 is None and mol2 is None:
        logging.warning(f"None mol returned for smiles: {smi1}, {smi2}")
        return None, None
    elif mol1 is None:
        logging.warning(f"None mol returned for smiles: {smi1}")
        return None, None
    elif mol2 is None:
        logging.warning(f"None mol returned for smiles: {smi2}")
        return None, None        

    res = rdFMCS.FindMCS([mol1, mol2], 
                            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            matchValences=do_valence, 
                            matchChiralTag=False,
                            ringMatchesRingOnly=True,
                            completeRingsOnly=True,
                            timeout=1
            )

    if res.canceled:
        logging.warning(f"MCS timeout for smiles: {smi1}, {smi2}")
        return None, None
    elif norm == 'min':
        mcs_idx, tot = res.numAtoms, min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
    elif norm == 'max':
        mcs_idx, tot = res.numAtoms, max(mol1.GetNumAtoms(), mol2.GetNumAtoms())
    
    return mcs_idx, tot  

def mcs_align(smarts1, smarts2, norm='max', do_valence=True):
    r1, p1 = smarts_to_sub_smiles(smarts1)
    r2, p2 = smarts_to_sub_smiles(smarts2)
    r1, r2, p1, p2 = [list(set(elt)) for elt in [r1, r2, p1, p2]] # Remove stoich
    
    if len(r1 + p1) != len(r2 + p2):
        logging.warning(f"Unequal number of unique substrates for {smarts1}, {smarts2}")
        return 0
    
    subs1 = r1 + p1
    rxn_mcs_idxs = []
    for rperm in permutations(r2):
        for pperm in permutations(p2):
            subs2 = rperm + pperm
            running_mcs, running_tot = 0, 0
            for elt in zip(subs1, subs2):
                this_mcs, this_tot = mcs_index(elt[0], elt[1], do_valence, norm)
                
                # Timeout or san issue
                if this_mcs is None and this_tot is None:
                    return 0
                
                running_mcs += this_mcs
                running_tot += this_tot

            if running_tot == 0:
                atom_ave_mcs = 0
            else:
                atom_ave_mcs = running_mcs / running_tot
            rxn_mcs_idxs.append(atom_ave_mcs)

    rxn_mcs = max(rxn_mcs_idxs)
    return rxn_mcs
            
def is_mcs_similar(x, putative_smarts, jni_up, rhead_id_to_smarts, threshold=0.9):
        rheas = [int(elt.split(":")[1]) for elt in jni_up[x]['rhea_ids']]
        bool_list = []
        for rid in rheas:
            above_threshold = False
            for rhea_smarts in rhea_id_to_smarts.get(rid, []):
                mcs_rxn_idx = mcs_align(putative_smarts, rhea_smarts)
                if mcs_rxn_idx > threshold:
                    above_threshold = True
                    break

            bool_list.append(above_threshold)

        return any(bool_list)

def filter_known_rxns_uniprot(filter_fcn, known_rxns, jni_up=jni_up,
                                rhea_id_to_rhashes=rhea_id_to_rhashes,
                                rhea_id_to_stoichless_rhashes=rhea_id_to_stoichless_rhashes,
                                rhea_id_to_smarts=rhea_id_to_smarts):
    '''
    Wraparound filter dict to process known rxns entries    
    '''
    temp = deepcopy(known_rxns)
    filtered_known_rxns = {}
    for k,v in tqdm(temp.items(), desc='known rxns'):
        if filter_fcn is is_rhash_equivalent:
            _filter_fcn = lambda x: filter_fcn(x, k, jni_up, rhea_id_to_rhashes) # Pass putative rhash, k and dicts

        elif filter_fcn is is_stoichless_rhash_equivalent:
            _filter_fcn = lambda x: filter_fcn(x, v['stoichless_rhash'], jni_up, rhea_id_to_stoichless_rhashes) # Pass stoichless putative rhash

        elif filter_fcn is is_mcs_similar:
            _filter_fcn = lambda x: filter_fcn(x, v['smarts'], jni_up, rhea_id_to_smarts) # Pass smarts, rhea smarts lookup

        else:
            _filter_fcn = filter_fcn

        ds_uids = list(filter(_filter_fcn, v['uniprot_ids']))
        filtered_known_rxns[k] = v
        filtered_known_rxns[k]['uniprot_ids'] = ds_uids
    return filtered_known_rxns

# Filter known reactions' uniprot ids
known_rxns_found_uid = filter_known_rxns_uniprot(in_uniprot, known_rxns)
n_any, n_rxns = ct_any_ids(known_rxns_found_uid)
print(f"{n_any} / {n_rxns} rxns with uniprot ids I found in uniprot")

known_rxns_found_rhea_rxn = filter_known_rxns_uniprot(has_rxn, known_rxns_found_uid)
n_any, n_rxns = ct_any_ids(known_rxns_found_rhea_rxn)
print(f"{n_any} / {n_rxns} rxns with uniprot with reactions")

known_rxns_w_rhash_eq = filter_known_rxns_uniprot(is_rhash_equivalent, known_rxns_found_rhea_rxn)
n_any, n_rxns = ct_any_ids(known_rxns_w_rhash_eq)
print(f"{n_any} / {n_rxns} rxns with uniprot with rhash equivalent rhea reactions")

unmatched_by_rhash = {}
for k in known_rxns.keys():
    if known_rxns_found_rhea_rxn[k]['uniprot_ids'] and not known_rxns_w_rhash_eq[k]['uniprot_ids']:
        unmatched_by_rhash[k] = known_rxns_found_rhea_rxn[k]

side_by_side = defaultdict(list)
rhea_rhashes_to_smarts = defaultdict(set)
for key, kr in unmatched_by_rhash.items():
    for uid in kr['uniprot_ids']:
        for rhea in jni_up[uid]['rhea_ids']:
            rhea = int(rhea.split(":")[1])
            if rhea in rhea_id_to_rhashes:
                rhea_hash = rhea_id_to_rhashes[rhea]
                rhea_smarts = rhea_id_to_smarts[rhea]
                side_by_side[(key, rhea_hash)].append(rhea)
                rhea_rhashes_to_smarts[rhea_hash].add(rhea_smarts[0])
                rhea_rhashes_to_smarts[rhea_hash].add(rhea_smarts[1])

matches_by_mcs = {} # Indexed in same way as side_by_side
for i, k in tqdm(enumerate(list(side_by_side.keys())[:200])):
    kr_rhash = k[0]
    rhea_rhash = k[1]
    kr_smarts = unmatched_by_rhash[kr_rhash]['smarts']
    
    mcses = []
    for sma in rhea_rhashes_to_smarts[rhea_rhash]:
        rxn_mcs = mcs_align(kr_smarts, sma, norm='max', do_valence=True)
        mcses.append(rxn_mcs)
    
    matches_by_mcs[(kr_rhash, rhea_rhash)] = max(mcses)

    if i % 500 == 0:
        print(f"Saving {i}th reaction-pair")
        with open("../data/mapping/mcs_similarity_known_rhea_rxn_pairs.pkl", 'wb') as f:
            pickle.dump(matches_by_mcs, f)

        

