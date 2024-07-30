from argparse import ArgumentParser
from src.utils import load_json
from src.operator_mapping import match_template, map_rxn2rule
from src.cheminfo_utils import standardize_smarts_rxn
import pandas as pd
import multiprocessing as mp
from itertools import chain

# Run from cmd
# E.g., python map_rxns.py JN3604IMT_rules.tsv swissprot_unmapped.json swissprot_unmapped.json

'''
Args
-----
rules - Path to operators tsv file w/ columns: Name | Reactants | SMARTS | Products
reactions - Path to reactions json w/ {unique_id: SMARTS} where SMARTS:str like 'reactant.reactant>>product.product'
output - Path to save mapping results

Returns
-------
Mapping results in a tsv w/ columns: Reaction ID | Rule | Aligned smarts | Reaction center
and every row is one Reaction-Rule map pair
'''
parser = ArgumentParser()
parser.add_argument("rules", help='Path to operators tsv file w/ columns: Name | Reactants | SMARTS | Products')
parser.add_argument("reactions", help="Path to reactions json w/ {unique_id: SMARTS} where SMARTS:str like 'reactant.reactant>>product.product'")
parser.add_argument("output", help="Path to save mapping results")
args = parser.parse_args()
do_template = True # Whether to enforce template matching, ie cofactors
return_rc = True # Whether to return reaction center while mapping operators
pre_standardized = True
rm_stereo = True

def map_reaction(
        rxn_id,
        sansmarts,
        rules,
        smi2paired_cof,
        smi2unpaired_cof,
        return_rc
        ):
    output_rows = []
    for rule_name, row in rules.iterrows(): # Iterate over rules
        did_map = False
        rule_reactants_template, rule_smarts, rule_products_template = row

        matched_idxs = match_template(sansmarts, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof)

        if len(matched_idxs) > 0:
            res = map_rxn2rule(sansmarts, rule_smarts, return_rc=return_rc, matched_idxs=matched_idxs) # Map if have template matches
            did_map = res['did_map']

        if did_map:
            print(f"{rxn_id} => {rule_name}")
            output_rows.append([rxn_id, rule_name, res['aligned_smarts'], res['reaction_center']])
    
    return output_rows

def try_standardize_smarts_rxn(smarts, rm_stereo):
    try:
        return standardize_smarts_rxn(smarts, remove_stereo=rm_stereo)
    except:
        return None
    
def prep_reaction(
        rxn_id,
        entry,
        pre_standardized,
        rm_stereo
    ):
    smarts = entry['smarts']
    if not pre_standardized:
        sansmarts = try_standardize_smarts_rxn(smarts, rm_stereo)
        
        if sansmarts is None:
            print(f"Unable to sanitize reaction: {rxn_id} w/ SMARTS: {smarts}")
        
        return sansmarts
    else:
        return smarts

def process_reaction(
        pair,
        rules,
        smi2paired_cof,
        smi2unpaired_cof,
        return_rc,
        pre_standardized,
        rm_stereo,
    ):
    rxn_id, entry = pair
    sansmarts = prep_reaction(
        rxn_id=rxn_id,
        entry=entry,
        pre_standardized=pre_standardized,
        rm_stereo=rm_stereo
    )

    output_rows = map_reaction(
        rxn_id=rxn_id,
        sansmarts=sansmarts,
        rules=rules,
        smi2paired_cof=smi2paired_cof,
        smi2unpaired_cof=smi2unpaired_cof,
        return_rc=return_rc
    )
    return output_rows

# Read in rules
rules = pd.read_csv(args.rules, sep='\t')
rules.set_index("Name", inplace=True)
rules.drop('Comments', axis=1, inplace=True)

rxns = load_json(args.reactions) # Read in reactions
n_rxns = len(list(rxns.keys())) # Total no. reactions to map

# Read in cofactor lookup tables
smi2unpaired_cof = load_json('/home/stef/bottle/data/cofactors/smi2unpaired_cof.json')
smi2paired_cof = load_json('/home/stef/bottle/data/cofactors/smi2paired_cof.json')
smi2paired_cof = {tuple(k.split(",")): tuple(v.split(",")) for k,v in smi2paired_cof.items()}

output_cols = ["Reaction ID", "Rule", "Aligned smarts", "Reaction center"]
output_data = []
processed_reactions = []
mapped_rules = []
reaction_centers = []
pairs = list(rxns.items())

def process_pair(pair):
    output_rows = process_reaction(
        pair,
        rules,
        smi2paired_cof,
        smi2unpaired_cof,
        return_rc,
        pre_standardized,
        rm_stereo
    )
    return output_rows

with mp.Pool() as pool:
    res = pool.map(process_pair, pairs)
    
    
output_data = list(chain(*res)) 
   

df = pd.DataFrame(
    data=output_data,
    columns=output_cols
)

df.to_csv(
    path_or_buf=args.output,
    sep='\t',
    index=False
)