from src.utils import sort_x_by_y
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rxnmapper import RXNMapper
rxnmapper = RXNMapper()

def atom_map(rxn_sma):
    return rxnmapper.get_attention_guided_atom_maps([rxn_sma])[0]['mapped_rxn']

def get_sub_mol(mol, atom_idxs):
    '''
    Given mol object and tuple
    of atom indices returns mol 
    of the part of the mol with
    atoms at those indices incl bonds if applicable
    '''
    sub_mol = Chem.RWMol() # New mol object to be edited
    old2new = {} # {idx_in_mol:idx_in_sub_mol}
    
    # Add atoms, set formal charge
    for i, elt in enumerate(atom_idxs):
        atom = mol.GetAtomWithIdx(elt)
        fc = atom.GetFormalCharge()
        amap_num = atom.GetAtomMapNum()
        anum = atom.GetAtomicNum()
        old2new[elt] = i
        sub_mol.AddAtom(Chem.Atom(anum))
        sub_mol.GetAtomWithIdx(i).SetFormalCharge(fc)
        sub_mol.GetAtomWithIdx(i).SetAtomMapNum(amap_num)

    # Add bonds
    for a1 in atom_idxs:
        for a2 in atom_idxs:
            if a1 < a2:
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None:
                    b_type = bond.GetBondType() # Single, double, ...
                    sub_mol.AddBond(old2new[a1], old2new[a2], b_type)

    # Trust in Greg Landrum
    sub_mol = Chem.Mol(sub_mol)
    Chem.SanitizeMol(sub_mol)

    return sub_mol

def align_substrates(rcs, remaining):
    '''
    Aligns substrates of two reactions based on 
    their reactions centers being equivalent.
    Args
        - rcs: List of lists of mol obj rxn ctrs
        - remaining: List of lists of remaining rxn ctr idxs
    Returns
        - List of aligned idx pairs (rxn1_idx, rxn2_idx)
        or None if 
    '''

    for i in remaining[0]:
        for j in remaining[1]:
            rc1 = rcs[0][i]
            rc2 = rcs[1][j]
            res = rdFMCS.FindMCS([rc1, rc2], bondCompare=rdFMCS.BondCompare.CompareOrder,
               atomCompare=rdFMCS.AtomCompare.CompareElements, matchValences=True
               )
            if (res.numAtoms == max(rc1.GetNumAtoms(), rc2.GetNumAtoms())) &\
                (res.numBonds == max(rc1.GetNumBonds(), rc2.GetNumBonds())):
                idx_pair = (i, j)
                return idx_pair
            
    return None

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

def align_atom_map_nums(rxns, rcs, rc_atoms):
    '''
    Re-label atom map #'s in substrate pairs
    so that reaction centers of aligned subs
    take atom map #'s from rxn1 subs (arbitrarily).
    This is necessary for PRC MCS. Substrates should
    already be aligned by virtue of having matching
    rxn ctrs.
    Args
        - rxns: List of rxn objects [rxn1, rxn2]
        - rcs: List of lists of rxn ctr mol objects
        for the substrates of rxn1 & rxn2
        - rc_atoms: List of tuples of reacting atom
        idxs
    Returns
        - rxns: New list of rxn objects with re-labeled
        atom map numbers for rxn2 substrates
    '''
    rxns = rxns.copy() # Defensive copy; avoid side effects
    for j, rc1 in enumerate(rcs[0]): # For rxn ctr from rxn 1
        mol2 = rxns[1].GetReactantTemplate(j) # Get corresponding substrate from rxn2
        ratoms2 = rc_atoms[1][j] # And rxn ctr atom idxs of this molecule
        rc1_atom_idxs = [i for i in range(len(rc1.GetAtoms()))]
        
        # Get atom map #'s from rxn1 rxn ctr
        rc_amap_nums = []
        for elt in rc1_atom_idxs:
            rc_amap_nums.append(rc1.GetAtomWithIdx(elt).GetAtomMapNum())

        # Property hashes of rxn ctr from rc1 and mol2
        # give true id of each atom
        mol2_prop_hashes = get_property_hashes(mol2, ratoms2)
        rc1_prop_hashes = get_property_hashes(rc1, rc1_atom_idxs)

        # Sort rxn ctr idxs from rxn2 sub and atom map #'s from 
        # rxn1 sub by their true identity
        ratoms2, mol2_prop_hashes = tuple(sort_x_by_y(ratoms2, mol2_prop_hashes))
        rc_amap_nums, rc1_prop_hashes = tuple(sort_x_by_y(rc_amap_nums, rc1_prop_hashes))
        
        assert mol2_prop_hashes == rc1_prop_hashes

        # Re-label mol from rxn2
        for i, elt in enumerate(ratoms2):
            mol2.GetAtomWithIdx(elt).SetAtomMapNum(rc_amap_nums[i])

    return rxns

def get_prc_mcs(rxns, rcs, rc_atoms, norm='min atoms'):
    '''
    Returns peri-reaction-center maximum common
    substructure for each aligned pair of substrates
    Args
        - rxns: List of rxn objects [rxn1, rxn2]
        - rcs: List of lists of rxn ctr mol objects
        for the substrates of rxn1 & rxn2
        - rc_atoms: List of tuples of reacting atom
        idxs
        - norm: Normalization to get an index out of
        prcmcs. 'min atoms' (default) normalizes by smaller
        of the two substrates, 'max atoms' by the larger
    Returns
        - prc_mcs: List of prc_mcs index values
    '''
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
    for rc in rcs[0]:
        for atom in rc.GetAtoms():
            atom.SetIsotope(atom.GetAtomMapNum() * atom.GetAtomicNum() * 99)

    # Get prc mcs
    prc_mcs = []
    for i in range(rxns[0].GetNumReactantTemplates()):
        subs = [rxns[0].GetReactantTemplate(i), rxns[1].GetReactantTemplate(i)]
        rc = rcs[0][i] # Rxn ctr from rxn1 has right atom map # for FindMCS seed

        res = rdFMCS.FindMCS(subs, seedSmarts=Chem.MolToSmarts(rc),
                                atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                                bondCompare=rdFMCS.BondCompare.CompareOrder
                            ) # Left out matchValences for now... unk issue
        
        # Compute prc mcs index
        if norm == 'min atoms':
            prc_mcs.append(res.numAtoms / min(subs[0].GetNumAtoms(), subs[1].GetNumAtoms()))
        elif norm == 'max atoms':
            prc_mcs.append(res.numAtoms / max(subs[0].GetNumAtoms(), subs[1].GetNumAtoms()))
    
    return prc_mcs


# Tests
if __name__ == '__main__':
    # Test that align atom map nums
    pass

