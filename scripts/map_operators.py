from argparse import ArgumentParser
from src.utils import load_json
from src.operator_mapping import match_template, map_rxn2rule, expand_paired_cofactors, expand_unpaired_cofactors
from src.cheminfo_utils import standardize_smarts_rxn
import pandas as pd
import multiprocessing as mp
from itertools import chain
from tqdm import tqdm

def try_standardize_smarts_rxn(smarts):
    try:
        return standardize_smarts_rxn(smarts)
    except:
        return None

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

def prep_reaction(
        rxn_id,
        entry,
        pre_standardized,
        rm_stereo
    ):
    smarts = entry['smarts']
    if pre_standardized:
        return smarts
    else:
        sansmarts = try_standardize_smarts_rxn(smarts, rm_stereo)
        
        if sansmarts is None:
            print(f"Unable to sanitize reaction: {rxn_id} w/ SMARTS: {smarts}")
        
        return sansmarts

def process_reaction(
        rxn_id,
        rxn_entry,
        rules,
        smi2paired_cof,
        smi2unpaired_cof,
        return_rc,
        pre_standardized,
        rm_stereo,
    ):
    sansmarts = prep_reaction(
        rxn_id=rxn_id,
        entry=rxn_entry,
        pre_standardized=pre_standardized,
        rm_stereo=rm_stereo
    )

    if sansmarts is None: # Unable to standardize
        return []

    output_rows = map_reaction(
        rxn_id=rxn_id,
        sansmarts=sansmarts,
        rules=rules,
        smi2paired_cof=smi2paired_cof,
        smi2unpaired_cof=smi2unpaired_cof,
        return_rc=return_rc
    )
    return output_rows

def mp_wrap(bunch):
    return process_reaction(**bunch)

if __name__ == '__main__':
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
    return_rc = False # Whether to return reaction center while mapping operators
    pre_standardized = True
    rm_stereo = True
    k_tautomers = 10

    # Read in rules
    rules = pd.read_csv(args.rules, sep='\t')
    rules.set_index("Name", inplace=True)
    rules.drop('Comments', axis=1, inplace=True)

    rxns = load_json(args.reactions) # Read in reactions
    n_rxns = len(list(rxns.keys())) # Total no. reactions to map

    # Read in cofactor lookup tables
    paired_ref = pd.read_csv('../data/cofactors/paired_cofactors_reference.tsv', sep='\t')
    unpaired_ref = pd.read_csv('../data/cofactors/unpaired_cofactors_reference.tsv', sep='\t')
    smi2paired_cof = expand_paired_cofactors(paired_ref, k=k_tautomers)
    smi2unpaired_cof = expand_unpaired_cofactors(unpaired_ref, k=k_tautomers)

    output_cols = ["Reaction ID", "Rule", "Aligned smarts", "Reaction center"]
    output_data = []
    processed_reactions = []
    mapped_rules = []
    reaction_centers = []
    bunches = [
        {
            'rxn_id':k,
            'rxn_entry':v,
            'rules':rules,
            'smi2paired_cof':smi2paired_cof,
            'smi2unpaired_cof':smi2unpaired_cof,
            'return_rc':return_rc,
            'pre_standardized':pre_standardized,
            'rm_stereo':rm_stereo
        }
        for k,v in rxns.items()
    ]

    with mp.Pool() as pool:
        res = list(tqdm(pool.imap(mp_wrap, bunches), total=len(bunches)))

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