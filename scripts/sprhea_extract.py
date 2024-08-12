import pandas as pd
from src.utils import load_json, save_json
from src.post_processing import Enzyme
from src.cheminfo_utils import standardize_smiles, clean_up_rhea_rxn
from math import isnan
from collections import defaultdict
from itertools import chain

def sanitize_smiles(smiles):
    '''
    Wrap standardize smiles in try except
    '''
    try:
        return standardize_smiles(smiles)
    except:
        return None
    
def extract_enzyme(upid, up_pull):
    '''
    Pull enzyme info from dataframe up_pull at 
    index upid
    '''
    row = up_pull.loc[upid]
    keys = [
        'uniprot_id',
        'sequence',
        'existence',
        'reviewed',
        'ec',
        'organism',
        'name'
    ]
    cols = [
        'Sequence',
        'Protein existence',
        'Reviewed',
        'EC number',
        'Organism',
        'Protein names'
    ]
    args = [upid]
    for col in cols:
        if type(row[col]) is float:
            if isnan(row[col]):
                args.append(None)
            else:
                args.append(row[col])
        else:
            args.append(row[col])

    e = Enzyme(**dict(zip(keys, args)))
    return e.to_dict()
    
def extract_smi2name(rheas, sansmarts, rhea2smi2name):
    '''
    Construct smiles 2 names look ups for each smiles in 
    sansmarts. Pulls info from rhea2smi2name. Appends names
    when there are multiple rhea ids for one san smarts
    '''
    smiles = list(chain(*[elt.split('.') for elt in sansmarts.split('>>')]))
    combosmi2name = {}
    for rhea in rheas:
        for smi, name in rhea2smi2name[rhea].items():
            if smi in combosmi2name:
                combosmi2name[smi] += f";{name}" # Append multiple names
            else:
                combosmi2name[smi] = name
    
    # Remove duplicate names
    for k,v in combosmi2name.items():
        simple_v = ";".join(set(v.split(';')))
        combosmi2name[k] = simple_v
    
    return {smi: combosmi2name.get(smi, None) for smi in smiles}
         

'''
Load in:
1. Rhea-SMARTS
2. UniProt pull
3. All Rhea-UniProt pairs found in pull and in Rhea flat files
4. Rhea to SMILES2Names lookup
'''
rhea_smarts = pd.read_csv('/home/stef/bottle/data/sprhea/rhea-reaction-smiles.tsv', '\t', header=None)
rhea_smarts.columns = ["rhea_id", "smarts"]
rhea_smarts.set_index('rhea_id', inplace=True)

up_pull = pd.read_csv('/home/stef/bottle/data/sprhea/uniprotkb_reviewed_true_AND_proteins_wi_2024_02_29.tsv', sep='\t')
up_pull.set_index("Entry", inplace=True)

rhea2upid = load_json("/home/stef/bottle/data/sprhea/rhea2uniprot_all.json")

rhea_smiles_names = load_json("/home/stef/bottle/data/sprhea/rhea_smiles_names.json")

rhea2upid = {int(k): v for k,v in rhea2upid.items()}
rhea_smiles_names = {int(k): v for k,v in rhea_smiles_names.items()}
all_rheas = [rid for rid in rhea2upid.keys()]

# Construct Rhea ID to SMARTS dict
rhea2smarts = {rid: rhea_smarts.loc[rid, 'smarts'] for rid in all_rheas}

# Construct Rhea ID to SMILES to name dict
rhea2smi2name = defaultdict(lambda : defaultdict())
for rid, elt in rhea_smiles_names.items():
    if elt is None:
        continue
    for side in ['left', 'right']:
        for pair in elt[side]:
            rhea2smi2name[rid][pair['smiles']] = pair['label']

# Initialize SMILES to sanitized SMILES dict
smi2sansmi = {}

# Sanitize Rhea ID to SMILES
print("Sanitizing smiles in smiles2name lookup")
rhea2sansmi2name = defaultdict(lambda : defaultdict())
for rid, elt in rhea2smi2name.items():
    for smi, name in elt.items():
        if smi in smi2sansmi:
            rhea2sansmi2name[rid][smi2sansmi[smi]] = name
        else:
            sansmi = sanitize_smiles(smi)
            rhea2sansmi2name[rid][sansmi] = name
            smi2sansmi[smi] = sansmi

# Iterate over Rhea IDs w/ SMARTS 
# Sanitize SMARTS SMILES by SMILES, pulling from SMILES to san SMILES where possible. 
# Store Rhea ID to {san SMARTS, smi2name}.
# Store rhash to Rhea ID.
# Store san SMARTS to rxnid. #
print("Sanitizing Smarts")
rxnid2rhea = defaultdict(list)
sansmarts2rxnid = {}
rxn_no = 0
for rid, smarts in rhea2smarts.items():
    cu_smarts = clean_up_rhea_rxn(smarts) # MUST remove numbered asterisks before sanitizing smiles
    reactants, products = [elt.split('.') for elt in cu_smarts.split('>>')]
    san_reactants, san_products = [], []
    for r in reactants:
        if r in smi2sansmi:
            san_reactants.append(smi2sansmi[r])
        else:
            sanr = sanitize_smiles(r)
            san_reactants.append(sanr)
            smi2sansmi[r] = sanr
    for p in products:
        if p in smi2sansmi:
            san_products.append(smi2sansmi[p])
        else:
            sanp = sanitize_smiles(p)
            san_products.append(sanp)
            smi2sansmi[p] = sanp

    if None in san_reactants or None in san_products: # Skip if standardize_smiles fails
        continue
    
    sansmarts = ".".join(sorted(san_reactants)) + '>>' + ".".join(sorted(san_products))

    if sansmarts in sansmarts2rxnid:
        rxnid2rhea[sansmarts2rxnid[sansmarts]].append(rid)
    else:
        sansmarts2rxnid[sansmarts] = rxn_no
        rxnid2rhea[rxn_no].append(rid)
        rxn_no += 1

# Pull san SMARTS, Enzymes.to_dict(), Rhea IDs, smi2name, reverse fields
# into rxn-indexed dict #
known_rxns = {}
for sansmarts, rxnid in sansmarts2rxnid.items():
    rev_smarts = sansmarts.split('>>')[1] + ">>" + sansmarts.split('>>')[0]
    rev_rxnid = sansmarts2rxnid.get(rev_smarts, None)

    this_rheas = [rhea for rhea in rxnid2rhea[rxnid]]
    
    enzymes = []
    upids = set(chain(*[rhea2upid[elt] for elt in this_rheas]))
    for upid in upids:
        if upid not in up_pull.index:
            continue
        else:
            enzymes.append(extract_enzyme(upid, up_pull))

    smi2name = extract_smi2name(this_rheas, sansmarts, rhea2sansmi2name)

    known_rxns[rxnid] = {
        'smarts':sansmarts,
        'min_rule':None,
        'imt_rules':None,
        'smi2name':smi2name,
        'enzymes':enzymes,
        'reaction_center':None,
        'reverse':rev_rxnid,
        'rhea_ids':this_rheas,

    }

# Remove transport reactions
trans = []
for k,v in known_rxns.items():
    if v['reverse'] == k:
        trans.append(k)

for k in trans:
    known_rxns.pop(k, None)

save_json(known_rxns, "/home/stef/bottle/data/sprhea/sprhea_240310_v3.json")