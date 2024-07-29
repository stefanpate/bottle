import json
import numpy as np
from itertools import permutations, product, chain
from rdkit import Chem
from rdkit.Chem import AllChem, CanonSmiles
import re

def map_rxn2rule(rxn, rule, return_rc=False, matched_idxs=None, max_products=10000):
    '''
    Maps reactions to SMARTS-encoded reaction rule.
    Args:
        - rxn: List of lists, each sublist with smiles of 
        substrates with the correct multiplicity / stoichiometry
        - rule: smarts string
        - return_rc: Return reaction center
        - matched_idxs: Indices of reaction reactants in the order they match the smarts
        reactants templates
    Returns:
        - did_map (bool)
        - rc: Tuple of w/ reaction center atoms for each reactant or empty tuple
    '''
    did_map = False
    reactants, products = rxn
    
    # If neither missing nor unparseable smiles, continue with mapping
    products = sorted(products)
    operator = Chem.rdChemReactions.ReactionFromSmarts(rule) # Make reaction object from smarts string
    reactants_mol = [Chem.MolFromSmiles(elt) for elt in reactants] # Convert reactant smiles to mol obj
    rule_substrate_cts = count_reactants(rule) # [n_reactants, n_products] in a rule
    rxn_substrate_cts = [len(reactants), len(products)]

    # Check if number of reactants / products strictly match
    # rule to reaction. If not return false
    if rule_substrate_cts != rxn_substrate_cts:
        return did_map, tuple()
    
    # If not enforcing templates,
    # get all permutations of reactant
    # indices
    if matched_idxs is None:
        matched_idxs = list(permutations([i for i in range(len(reactants))]))
        
    # For every permutation of that subset of reactants
    for idx_perm in matched_idxs:
        perm = tuple([reactants_mol[idx] for idx in idx_perm]) # Re-order reactants based on allowable idx perms
        outputs = operator.RunReactants(perm, maxProducts=max_products) # Apply rule to that permutation of reactants

        did_map = compare_operator_outputs_w_products(outputs, products)

        if did_map:
            break # out of permutations-of-matched-idxs loop

    if did_map and not return_rc: # Mapped and don't want rc
        return did_map, tuple()

    elif did_map and return_rc: # Mapped and want rc
        patts = get_lhs_patts_from_operator(rule)
        patts = [Chem.MolFromSmarts(elt) for elt in patts]

        if len(patts) != len(perm):
            raise Exception("Something wrong. There should be same number of operator fragments as reaction reactants") # TODO This the right way to raise exceptions?
        
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
                
                # Re-order rcs back to order of original reaction smiles
                rc = [None for elt in smc]
                for i, idx in enumerate(idx_perm):
                    rc[idx] = smc[i]

                return did_map, tuple(rc)
            
            # Deprotect & try again
            for j, reactant_rc in enumerate(smc):
                all_but = all_putative_rc_atoms[j] - set(reactant_rc) # To protect: "all but current rc"
                for protect_idx in all_but:
                    perm[j].GetAtomWithIdx(protect_idx).ClearProp('_protected')

    else: # Did not map
        return did_map, tuple()

def match_template(rxn, rule_reactants_template, rule_products_template, smi2paired_cof, smi2unpaired_cof):
    '''
    Returns the permuted indices corresponding to
    a match between reactant and rule templates
    '''
    reactants_smi, products_smi = rxn
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

def count_reactants(rule_smarts):
    '''
    Counts number of reactants in a SMARTS-
    encoded operator
    '''
    sides = rule_smarts.split('>>')
    cts = []
    for side in sides:
        dot_split = side.split('.') # Reactants separated by '.'

        # But must catch where pieces of a single compound
        # are split by '.', in which case they'll be surrounded by ()
        left_split_parens = []
        right_split_parens = []
        for i, elt in enumerate(dot_split):
            if (elt[0] == '(') & (elt[-1] != ')'):
                left_split_parens.append(i)
            elif (elt[0] != '(') & (elt[-1] == ')'):
                right_split_parens.append(i)
                
        left_split_parens, right_split_parens = np.array(left_split_parens), np.array(right_split_parens)
        overcount = (right_split_parens - left_split_parens).sum()

        n = len(dot_split) - overcount
        cts.append(int(n))
    return cts

def get_lhs_patts_from_operator(smarts_str):

    # lhs smarts pattern
    lhs_smarts = smarts_str.split('>>')[0]
    lhs_smarts = re.sub(r':[0-9]+]', ']', lhs_smarts)

    # identify each fragment
    smarts_list = []
    temp_fragment = []

    # append complete fragments only
    for fragment in lhs_smarts.split('.'):
        temp_fragment += [fragment]
        if '.'.join(temp_fragment).count('(') == '.'.join(temp_fragment).count(')'):
            smarts_list.append('.'.join(temp_fragment))
            temp_fragment = []

            # remove component grouping for substructure matching
            if '.' in smarts_list[-1]:
                smarts_list[-1] = smarts_list[-1].replace('(', '', 1)[::-1].replace(')', '', 1)[::-1]

    return smarts_list

def compare_operator_outputs_w_products(outputs, products):
    for output in outputs:
        try:
            output = [CanonSmiles(Chem.MolToSmiles(elt)) for elt in output] # Convert pred products to canonical smiles
        except:
            output = [Chem.MolToSmiles(elt) for elt in output]
        
        output = sorted(output)

        # Compare predicted to actual products. If mapped, return True
        if output == products: 
            return True
        
        # Last, try fixing kekulization issues
        postsan_output = postsanitize_smiles(output)
        for elt in postsan_output: # Iterate over sets of outputs w/ diff tautomers
            if sorted(elt) == products:
                return True
            
    return False
