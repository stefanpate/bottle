import json
from rdkit import Chem
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

def shuffle_mol(mol):
    '''
    '''
    idxs = [i for i in range(len(list(mol.GetAtoms())))]
    random.shuffle(idxs)
    shuffle_mol = Chem.RWMol() # New mol object to be edited
    old2new = {}
    # Add atoms, set formal charge
    for i, elt in enumerate(idxs):
        atom = mol.GetAtomWithIdx(elt)
        fc = atom.GetFormalCharge()
        amap_num = atom.GetAtomMapNum()
        anum = atom.GetAtomicNum()
        old2new[elt] = i
        shuffle_mol.AddAtom(Chem.Atom(anum))
        shuffle_mol.GetAtomWithIdx(i).SetFormalCharge(fc)
        shuffle_mol.GetAtomWithIdx(i).SetAtomMapNum(amap_num)

    # Add bonds
    for a1 in idxs:
        for a2 in idxs:
            if a1 < a2:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None:
                    b_type = bond.GetBondType() # Single, double, ...
                    shuffle_mol.AddBond(old2new[a1], old2new[a2], b_type)

    shuffle_mol = Chem.Mol(shuffle_mol)
    Chem.SanitizeMol(shuffle_mol) # Trust in Greg Landrum
    return shuffle_mol

def rm_atom_map_num(smarts):    
    rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)
    for elt in rxn.GetReactants():
        for atom in elt.GetAtoms():
            atom.SetAtomMapNum(0)
    for elt in rxn.GetProducts():
        for atom in elt.GetAtoms():
            atom.SetAtomMapNum(0)
    unmapped_smarts = AllChem.ReactionToSmiles(rxn)
    return unmapped_smarts