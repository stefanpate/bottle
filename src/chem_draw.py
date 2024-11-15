from rdkit import Chem
from rdkit.Chem import Draw
import svgutils.transform as st
from svgutils.compose import Unit
from collections import Counter
from src.config import filepaths
import numpy as np
from copy import copy

plus = st.fromfile(filepaths['artifacts'] / 'plus.svg').to_str().decode()
arrow = st.fromfile(filepaths['artifacts'] / 'arrow.svg').to_str().decode()

def _add_elt(img: str, movex: int, rxn_img: st.SVGFigure):
        elt = st.fromstring(img)
        root = elt.getroot()
        root.moveto(movex, 0)
        rxn_img.append(root)
        movex += int(elt.width.strip('px'))
        return rxn_img, movex

def _draw_side(side: dict, movex: int, rxn_img: st.SVGFigure, plus: str, auto_scl: bool):
    for i, (smi, stoich) in enumerate(side.items()):
        img = draw_molecule(smi, stoich, auto_scl=auto_scl)
        rxn_img, movex = _add_elt(img, movex, rxn_img)

        if i < len(side) - 1:
            rxn_img, movex = _add_elt(plus, movex, rxn_img)

    return rxn_img, movex


def draw_reaction(rxn_sma: str, plus: str = plus, arrow: str = arrow, auto_scl: bool = False) -> st.SVGFigure:
 
    reactants, products = [Counter(elt.split('.')) for elt in rxn_sma.split('>>')]

    movex = 0
    rxn_img = st.SVGFigure()
    
    rxn_img, movex = _draw_side(reactants, movex, rxn_img, plus, auto_scl)
    rxn_img, movex = _add_elt(arrow, movex, rxn_img)
    rxn_img, movex = _draw_side(products, movex, rxn_img, plus, auto_scl)
    
    elt = st.fromstring(arrow)
    height = float(elt.height.strip('px')) # Assumed same for all elements
    rxn_img.width = Unit(movex)
    rxn_img.height = Unit(height)

    return rxn_img

def draw_molecule(smiles: str, stoich : int = 1, size: tuple = (200, 200), hilite_atoms : tuple = tuple(), auto_scl: bool = False):
    '''
    Draw molecule.

    Args
    ----
    mol:str
        Molecule SMILES
    stoich:int
        Stoichiometric coefficient
    size:tuple
        (width, height)
     hilite_atoms:tuple
        Atom indices to highlight
    '''
    mol = Chem.MolFromSmiles(smiles)
    
    # Catch failed MolFromSmiles
    if mol is None: 
        mol = Chem.MolFromSmiles(smiles, sanitize=False)

    if auto_scl:
        na = mol.GetNumAtoms()
        width = int(np.log10(na) * 200) + 25
        size = (width, size[-1])

    drawer = Draw.MolDraw2DSVG(*size)

    if stoich == 1:
        drawer.DrawMolecule(mol, highlightAtoms=hilite_atoms)
    else:
        drawer.DrawMolecule(mol, legend=f"({stoich})", highlightAtoms=hilite_atoms)
    
    drawer.FinishDrawing()
    img = drawer.GetDrawingText()

    return img

if __name__ == '__main__':
    smi = 'CC(=O)O'
    mimg = draw_molecule(smi)

    rxn_sma = 'CCC.CCC.CCO>>CCCCCCCCO'
    fig = draw_reaction(rxn_sma, plus, arrow)