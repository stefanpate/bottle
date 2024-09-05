from src.cheminfo_utils import standardize_mol, tautomer_expand, standardize_smarts_rxn
from itertools import permutations, product, chain
from rdkit import Chem
import re
import pandas as pd


def post_standardize(mol, do_canon_taut):
    '''
    Post-reaction standardization
        - Skip neutralization because assume operators only affect heavy atoms, not hydrogens and therefore
        protonation states
        - Skip find parent because assume I am not producing salts / fragments (TODO: pressure test this
        assumption)
    -
    '''
    do_neutralize = False
    do_find_parent = False
    
    return Chem.MolToSmiles(standardize_mol(mol, do_canon_taut=do_canon_taut, do_neutralize=do_neutralize, do_find_parent=do_find_parent))

def split_reaction(rxn_smarts):
    return tuple([elt.split(".") for elt in rxn_smarts.split(">>")])

def map_rxn2rule(rxn, rule, return_rc=False, matched_idxs=None, max_products=10000):
    '''
    Maps reactions to SMARTS-encoded reaction rule.
    Args:
        - rxn: Reaction SMARTS string
        - rule: smarts string
        - return_rc: Return reaction center
        - matched_idxs: Indices of reaction reactants in the order they match the smarts
        reactants templates
    Returns:
        - res:dict{
            did_map:bool
            aligned_smarts:str | None
            reaction_center:Tuple[tuple] | None
        }
    '''
    res = {
        'did_map':False,
        'aligned_smarts':None,
        'reaction_center':None,
    }
    reactants, unsorted_products = split_reaction(rxn)
    
    products = sorted(unsorted_products) # Canonical ordering for later comparison
    operator = Chem.rdChemReactions.ReactionFromSmarts(rule) # Make reaction object from smarts string
    reactants_mol = [Chem.MolFromSmiles(elt) for elt in reactants] # Convert reactant smiles to mol obj
    rule_substrate_cts = [len(get_patts_from_operator_side(rule, i)) for i in range(2)] # [n_reactants, n_products] in a rule
    rxn_substrate_cts = [len(reactants), len(products)]

    # Check if number of reactants / products strictly match
    # rule to reaction. If not return false
    if rule_substrate_cts != rxn_substrate_cts:
        return res
    
    # If not enforcing templates,
    # get all permutations of reactant
    # indices
    if matched_idxs is None:
        matched_idxs = list(permutations([i for i in range(len(reactants))]))
        
    # For every permutation of that subset of reactants
    # TODO: What if there are multiple match idxs that product the right outputs?
    for idx_perm in matched_idxs:
        perm = tuple([reactants_mol[idx] for idx in idx_perm]) # Re-order reactants based on allowable idx perms
        outputs = operator.RunReactants(perm, maxProducts=max_products) # Apply rule to that permutation of reactants

        if compare_operator_outputs_w_products(outputs, products):
            res['did_map'] = True
            res['aligned_smarts'] = ".".join([reactants[idx] for idx in idx_perm]) + ">>" + ".".join(unsorted_products)
            break # out of permutations-of-matched-idxs loop

    if res['did_map'] and not return_rc: # Mapped and don't want rc
        return res

    elif res['did_map'] and return_rc: # Mapped and want rc
        patts = get_patts_from_operator_side(rule, 0)
        patts = [Chem.MolFromSmarts(elt) for elt in patts]

        if len(patts) != len(perm):
            raise Exception("Something wrong. There should be same number of operator fragments as reaction reactants")
        
        substruct_matches = [perm[i].GetSubstructMatches(patts[i]) for i in range(len(patts))]
        ss_match_combos = product(*substruct_matches) # All combos of putative rcs of n substrates
        all_putative_rc_atoms = [set(chain(*elt)) for elt in substruct_matches] # ith element has set of all putative rc atoms of ith reactant

        for smc in ss_match_combos:

            # Protect all but rc currently considered in each reactant
            for j, reactant_rc in enumerate(smc):
                all_but = all_putative_rc_atoms[j] - set(reactant_rc) # To protect: "all but current rc"
                for protect_idx in all_but:
                    perm[j].GetAtomWithIdx(protect_idx).SetProp('_protected', '1')

            outputs = operator.RunReactants(perm, maxProducts=max_products) # Run operator with protected atoms

            # If found match
            if compare_operator_outputs_w_products(outputs, products):
                res['reaction_center'] = smc
                return res
            
            # Deprotect & try again
            for j, reactant_rc in enumerate(smc):
                all_but = all_putative_rc_atoms[j] - set(reactant_rc) # To protect: "all but current rc"
                for protect_idx in all_but:
                    perm[j].GetAtomWithIdx(protect_idx).ClearProp('_protected')

    return res # Did not map or failed getting RC

def match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof):
    '''
    Returns the permuted indices corresponding to
    a match between reactant and rule templates
    '''
    reactants_smi, products_smi = split_reaction(rxn)
    rule_reactants_template = tuple(rule_reactants_template.split(';'))
    rule_products_template = tuple(rule_products_template.split(';'))
    matched_idxs = [] # Return empty if no matches found
    # First check the cardinality of reactants, products matches
    if (len(rule_reactants_template) == len(reactants_smi)) & (len(rule_products_template) == len(products_smi)):

        reactants_template = ['Any' for elt in reactants_smi]
        products_template = ['Any' for elt in products_smi]

        # Search for unpaired cofactors first
        for i, r in enumerate(reactants_smi):
            if r in smi2unpaired_cof:
                reactants_template[i] = smi2unpaired_cof[r]

        for i, p in enumerate(products_smi):
            if p in smi2unpaired_cof:
                products_template[i] = smi2unpaired_cof[p]

        # Search for paired cofactors
        # Only overwriting should be PPi/Pi as phosphate donor/acceptor
        for i, r in enumerate(reactants_smi):
            for j, p in enumerate(products_smi):
                if (r, p) in smi2paired_cof:
                    reactants_template[i] = smi2paired_cof[(r, p)][0]
                    products_template[j] = smi2paired_cof[(r, p)][1]
                elif (p, r) in smi2paired_cof:
                    reactants_template[i] = smi2paired_cof[(p, r)][1]
                    products_template[j] = smi2paired_cof[(p, r)][0]

        reactants_idx_template = [(elt, i) for i, elt in enumerate(reactants_template)]

        # First try to products templates
        product_template_match = False
        for perm in permutations(products_template):
            if perm == rule_products_template:
                product_template_match = True

        # If product templates match
        # find permutations of reactant template that match
        # rule template and keep the indices of those good permutations
        # Else return empty list
        if product_template_match:
            for perm in permutations(reactants_idx_template):
                this_template, this_idx = list(zip(*perm))
                if this_template == rule_reactants_template:
                    matched_idxs.append(this_idx)

    return matched_idxs

def get_patts_from_operator_side(smarts_str, side):

    # Side smarts pattern
    smarts = smarts_str.split('>>')[side]
    smarts = re.sub(r':[0-9]+]', ']', smarts)

    # identify each fragment
    smarts_list = []
    temp_fragment = []

    # append complete fragments only
    for fragment in smarts.split('.'):
        temp_fragment += [fragment]
        if '.'.join(temp_fragment).count('(') == '.'.join(temp_fragment).count(')'):
            smarts_list.append('.'.join(temp_fragment))
            temp_fragment = []

            # remove component grouping for substructure matching
            if '.' in smarts_list[-1]:
                smarts_list[-1] = smarts_list[-1].replace('(', '', 1)[::-1].replace(')', '', 1)[::-1]

    return smarts_list

def compare_operator_outputs_w_products(outputs, products):
    # Try WITHOUT tautomer canonicalization
    for output in outputs:
        try:
            output = sorted([post_standardize(mol, do_canon_taut=False) for mol in output]) # Standardize and sort SMILES
        except:
            continue

        # Compare predicted to actual products. If mapped, return True
        if output == products: 
            return True
        
    # Try WITH tautomer canonicalization
    try:
        products = sorted([post_standardize(mol, do_canon_taut=True) for mol in products])
    except:
        return False
    
    for output in outputs:
        try:
            output = sorted([post_standardize(mol, do_canon_taut=True) for mol in output]) # Standardize and sort SMILES
        except:
            continue

        # Compare predicted to actual products. If mapped, return True
        if output == products: 
            return True
            
    return False

def expand_paired_cofactors(df, k):
    smi2cof = {}
    for _, row in df.iterrows():
        smi_exp_1 = tautomer_expand(row["Smiles 1"], k)
        smi_exp_2 = tautomer_expand(row["Smiles 2"], k)
        for combo in product(smi_exp_1, smi_exp_2):
            smi2cof[combo] = (row["Class 1"], row["Class 2"])

    return smi2cof

def expand_unpaired_cofactors(df, k):
    smi2cof = {}
    for _, row in df.iterrows():
        smi_exp = tautomer_expand(row["Smiles"], k)
        for smi in smi_exp:
            smi2cof[smi] = row["Class"]

    return smi2cof

def standardize_template_map(
        rxn:str,
        rule_row:pd.Series,
        smi2paired_cof:dict,
        smi2unpaired_cof:dict,
        return_rc:bool,
        pre_standardized:bool,
    ):
    '''
    Convenience function to standardize, match cofactor templates, 
    and map a reaction to a rule.

    Args
    ----
    rxn:str
        Reaction smarts 'r1.r2>>p1.p2'
    rule_row:pd.Series
        Rule info with columns: SMARTS | Reactants | Products containing:
        reaction smarts | reactants cofactor template | products cofactor template
    smi2paired_cof:dict
        Paired cofactor lookup
    smi2unpaired_cor:dict
        Unpaired cofactor lookup
    return_rc:bool
        Return reaction center if True
    pre_standardize:bool
        Standardize reaction smarts before mapping if False

    Returns
    -------
    res['did_map']:bool
        Rule mapped reaction
    res['aligned_smarts']:str
        Reaction smarts w/ LHS aligned to operator
    res['reaction_center']:Tuple[Tuple]
        Tuple of tuple of reaction center atom indices
    '''
    if not pre_standardized:
        try:
            rxn = standardize_smarts_rxn(rxn)
        except:
            print(f"Unable to standardize reaction: {rxn}")
    
    rule = rule_row['SMARTS']
    rule_reactants_template = rule_row["Reactants"]
    rule_products_template = rule_row["Products"]
    matched_idxs = match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof)

    if len(matched_idxs) == 0:
        return False, None, None
    
    res = map_rxn2rule(rxn, rule, return_rc, matched_idxs)

    return res['did_map'], res['aligned_smarts'], res['reaction_center']