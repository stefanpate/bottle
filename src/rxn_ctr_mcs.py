from src.cheminfo_utils import sort_x_by_y
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from itertools import permutations, chain, product
import re
import numpy as np

def align_reactants(known_reaction_smarts, known_reaction_rc, min_rule):
    patts = get_lhs_patts_from_operator(min_rule)
    patts = [Chem.MolFromSmarts(elt) for elt in patts]
    known_reactants = [Chem.MolFromSmiles(elt) for elt in known_reaction_smarts.split('>>')[0].split('.')]
    known_reaction_rc = [tuple(elt) for elt in known_reaction_rc] # to compare w/ getsubstrucmatches output below

    rct_idxs = np.arange(len(known_reactants))
    for perm_idx in permutations(rct_idxs):
        perm = [patts[elt] for elt in perm_idx]
        matches = [mol.GetSubstructMatches(perm[i]) for i, mol in enumerate(known_reactants)]
        aligned = all([known_reaction_rc[i] in matches[i] for i in range(len(known_reaction_rc))])

        if aligned:
            return perm_idx


def get_property_hashes(mol, aidxs):
    '''
    Returns tuple of hashes of property vectors
    for atoms in a substructure of the mol object
    (indexed by aidxs). Property vectors include
    atomic #, formal charge, neighbor atoms (only
    those in the substructure), and bond types.
    '''

    prop_hashes = []
    for idx in aidxs:

        # Append atom info
        atom = mol.GetAtomWithIdx(idx)
        this_prop = [atom.GetAtomicNum(), atom.GetFormalCharge()]

        # Get bond & neighbor info
        bonds = []
        for bond in atom.GetBonds():

            # Neighbor is at "end atom idx" and part of the substructure
            if (bond.GetBeginAtomIdx() == idx) & (bond.GetEndAtomIdx() in aidxs):
                bonds.append((bond.GetBondType(), bond.GetEndAtom().GetAtomicNum()))
            
            # Neighbor is at "begin atom idx" and part of substructure
            if (bond.GetBeginAtomIdx() != idx) & (bond.GetBeginAtomIdx() in aidxs):
                bonds.append((bond.GetBondType(), bond.GetBeginAtom().GetAtomicNum()))

        if bonds:
            # Bonds have to be in consistent order
            bond_type, neighbor_atom_num = list(zip(*bonds))
            neighbor_atom_num, bond_type = sort_x_by_y(neighbor_atom_num, bond_type)
            bonds = list(zip(bond_type, neighbor_atom_num))

            this_prop += bonds # Append bond & neighbor info
        
        this_prop = tuple(this_prop)
        prop_hashes.append(hash(this_prop))

    prop_hashes = tuple(prop_hashes)

    return prop_hashes

def get_rc_mcs(
        rxns,
        rc_atoms,
        min_rule,
        norm='max atoms'
        ):
    '''
    Args
    ----
    rxns:List[str] - SMARTS of two reactant-aligned reactions 
    rc_atoms:List[Tuple[tuple]] - Idxs of reactant center atoms for each substrate
    min_rule:str - Min operator template shared btwn reactions
    norm:str - Normalization to get an index out of
        prcmcs. 'min atoms' normalizes by smaller
        of the two substrates, 'max atoms' by the larger
    Returns
    -------
    rc_mcs:List[float] - rc_mcs values per aligned substrate pair
    '''
    rxns = [AllChem.ReactionFromSmarts(sma, useSmiles=True) for sma in rxns]
    patts = get_lhs_patts_from_operator(min_rule)
    rc_fragments = [Chem.MolFromSmarts(patt) for patt in patts]

    # Set isotopes to customise FindMCS w/ atom map #
    for i, rxn in enumerate(rxns):
        for j, sub in enumerate(rxn.GetReactants()):
            for atom in sub.GetAtoms():
                atom_idx = atom.GetIdx()
                if atom_idx in rc_atoms[i][j]:
                    atom.SetIsotope(atom.GetAtomMapNum() * atom.GetAtomicNum() * 99) # Rxn ctr atom
                else:
                    atom.SetIsotope(atom.GetAtomicNum()) # Non rxn ctr atom

    # Set isotopes for 1st rxn's rxn ctrs same as rc within substrates
    for elt in rc_fragments:
        for atom in elt.GetAtoms():
            atom.SetIsotope(atom.GetAtomMapNum() * atom.GetAtomicNum() * 99)

    # Get prc mcs
    rc_mcs = []
    for i in range(rxns[0].GetNumReactantTemplates()):
        subs = [rxns[0].GetReactantTemplate(i), rxns[1].GetReactantTemplate(i)]

        for elt in subs:
            Chem.SanitizeMol(elt)
        
        frag = rc_fragments[i] # Rxn ctr from rxn1 has right atom map # for FindMCS seed

        res = rdFMCS.FindMCS(subs, seedSmarts=Chem.MolToSmarts(frag),
                                atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                                bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                matchChiralTag=False,
                                ringMatchesRingOnly=True,
                                completeRingsOnly=True,
                                matchValences=True,
                                timeout=10
                            )
        
        # Compute prc mcs index
        if res.canceled:
            rc_mcs.append(0)
        elif norm == 'min atoms':
            rc_mcs.append(res.numAtoms / min(subs[0].GetNumAtoms(), subs[1].GetNumAtoms()))
        elif norm == 'max atoms':
            rc_mcs.append(res.numAtoms / max(subs[0].GetNumAtoms(), subs[1].GetNumAtoms()))
    
    return tuple(rc_mcs)


def get_pred_rxn_ctr(pr_sma, min_rule):
    '''
    Returns atom idxs of LHS reaction center given predicted
    reaction smarts and the general operator template that generated it

    Args
    -----
    pr_sma:str
    rule:str
    (^^ both as 'sma1.sma2>>sma3.sma4')
    
    '''
    
    patts = get_lhs_patts_from_operator(min_rule)
    patts = [Chem.MolFromSmarts(elt) for elt in patts]
    reactants, products = [elt.split('.') for elt in pr_sma.split('>>')]
    rct_idxs = np.arange(len(reactants))
    for perm_idx in permutations(rct_idxs): # Note: Pickaxe doesn't keep PR reactants in order of operator template; must permute
        reactant_mols = [Chem.MolFromSmiles(reactants[idx]) for idx in perm_idx]
        operator = Chem.rdChemReactions.ReactionFromSmarts(min_rule) # Make reaction object from smarts string
        substruct_matches = [reactant_mols[i].GetSubstructMatches(patts[i]) for i in range(len(patts))]
        ss_match_combos = product(*substruct_matches) # All combos of putative rcs of n substrates
        all_putative_rc_atoms = [set(chain(*elt)) for elt in substruct_matches] # ith element has set of all putative rc atoms of ith reactant

        for smc in ss_match_combos:

            # Protect all but rc currently considered in each reactant
            for j, reactant_rc in enumerate(smc):
                all_but = all_putative_rc_atoms[j] - set(reactant_rc) # To protect: "all but current rc"
                for protect_idx in all_but:
                    reactant_mols[j].GetAtomWithIdx(protect_idx).SetProp('_protected', '1')

            outputs = operator.RunReactants(reactant_mols) # Run operator with protected atoms

            # If found match
            if compare_operator_outputs_w_products(outputs, products):
                rcts_perm = [reactants[idx] for idx in perm_idx]
                pr_sma_perm = ".".join(rcts_perm) + ">>" + ".".join(products)
                return smc, pr_sma_perm
            
            # Deprotect & try again
            for j, reactant_rc in enumerate(smc):
                all_but = all_putative_rc_atoms[j] - set(reactant_rc) # To protect: "all but current rc"
                for protect_idx in all_but:
                    reactant_mols[j].GetAtomWithIdx(protect_idx).ClearProp('_protected')
    
    return None, None

def compare_operator_outputs_w_products(outputs, products):
    products = sorted(products)
    for output in outputs:
        try:
            output = [Chem.CanonSmiles(Chem.MolToSmiles(elt)) for elt in output] # Convert pred products to canonical smiles
        except:
            output = [Chem.MolToSmiles(elt) for elt in output]
        
        output = sorted(output)

        # Compare predicted to actual products. If mapped, return True
        if output == products: 
            return True
            
    return False

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

def permute_by_idxs(to_permute, perm_idxs):
    '''
    Puts ith elt of to_permute in perm_idxs[i]th
    elt of perm. Returns perm
    '''
    perm = [None for _ in range(len(perm_idxs))]
    for i, pidx in enumerate(perm_idxs):
        perm[pidx] = to_permute[i]

    return perm

# Tests
if __name__ == '__main__':
    pass