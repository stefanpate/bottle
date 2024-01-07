import json
from rdkit import Chem
from rdkit.Chem import AllChem
import os

def ensure_dirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

def sort_x_by_y(x, y, reverse=False):
    yx = list(zip(y, x))
    sorted_yx = sorted(yx, reverse=reverse)
    y, x = list(zip(*sorted_yx))
    return x, y

def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)
    
def load_json(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data

def sanitize(list_of_smiles):
    '''
    Remove stereochem and canonicalize
    a list of smiles
    '''
    sanitized_smiles = []
    for elt in list_of_smiles:
        temp_mol = Chem.MolFromSmiles(elt)
        Chem.rdmolops.RemoveStereochemistry(temp_mol)
        sanitized_smiles.append(Chem.MolToSmiles(temp_mol))    
    return sanitized_smiles

def rxn_entry_to_smarts(rxn_entry):
    '''
    Convert our standard rxn json
    entry into a reaction smarts
    '''
    reactants = sanitize(list(rxn_entry[0].values()))
    products = sanitize(list(rxn_entry[1].values()))
    sma = ".".join(reactants) + ">>" + ".".join(products)
    return sma

def rm_atom_map_num(smarts):
    rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)

    # Remove atom map num and write mol smarts in order
    reactant_smas = []   
    for elt in rxn.GetReactants():
        for atom in elt.GetAtoms():
            atom.SetAtomMapNum(0)

        reactant_smas.append(Chem.MolToSmiles(elt))

    product_smas = []
    for elt in rxn.GetProducts():
        for atom in elt.GetAtoms():
            atom.SetAtomMapNum(0)

        product_smas.append(Chem.MolToSmiles(elt))

    unmapped_smarts = '.'.join(reactant_smas) + '>>' + '.'.join(product_smas)
    return unmapped_smarts