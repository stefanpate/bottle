from rdkit import Chem
from rdkit.Chem import Draw
import svgutils.transform as st
from svgutils.compose import Unit
from collections import Counter
import numpy as np
from pathlib import Path

root = Path(__file__).parent.parent

plus = st.fromfile(root / 'artifacts/plus.svg').to_str().decode()
arrow = st.fromfile(root / 'artifacts/arrow.svg').to_str().decode()
double_border = st.fromfile(root / 'artifacts/double_border.svg').to_str().decode()

def _add_elt(img: str, movex: int, rxn_img: st.SVGFigure):
        elt = st.fromstring(img)
        root = elt.getroot()
        root.moveto(movex, 0)
        rxn_img.append(root)
        movex += int(elt.width.strip('px'))
        return rxn_img, movex

def _draw_side(side: dict, movex: int, rxn_img: st.SVGFigure, plus: str, auto_scl: bool, sub_widths: list[int] = None):
    for i, (smi, stoich) in enumerate(side.items()):
        if sub_widths:
            img = draw_molecule(smi, stoich, auto_scl=auto_scl, size=(sub_widths[i], 200))
        else:
            img = draw_molecule(smi, stoich, auto_scl=auto_scl)
        
        rxn_img, movex = _add_elt(img, movex, rxn_img)

        if i < len(side) - 1:
            rxn_img, movex = _add_elt(plus, movex, rxn_img)

    return rxn_img, movex


def draw_reaction(rxn_sma: str, plus: str = plus, arrow: str = arrow, auto_scl: bool = False, size: tuple[int] = None) -> st.SVGFigure:
    '''
    Draw reaction

    Args
    ----
    rxn_sma:str
        Reaciton SMARTS R1.R2...>>P1.P2...
    plus:str
        SVG string of plus sign
    arrow:str
        SVG string of rightward reaction arrow
    auto_scl:bool
        If True, scales molecule image width proportional
        to log(# of atoms)

    Returns
    -------
    rxn_img:st.SVGFigure
        Which has save() and to_str() methods
    '''
 
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

def draw_molecule(molecule: str | Chem.Mol, stoich : int = 1, size: tuple = (200, 200), hilite_atoms : tuple = tuple(), auto_scl: bool = False) -> str:
    '''
    Draw molecule to svg string

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
    auto_scl:bool
        If True, scales molecule image width proportional
        to log(# of atoms)
    '''
    if isinstance(molecule, str):
        mol = Chem.MolFromSmiles(molecule)
    else:
        mol = molecule
    
    # Catch failed MolFromSmiles
    if mol is None: 
        mol = Chem.MolFromSmiles(molecule, sanitize=False)

    if auto_scl:
        na = mol.GetNumAtoms()
        width = int(np.log10(na) * size[0]) + 25
        size = (width, size[-1])

    drawer = Draw.MolDraw2DSVG(*[int(elt) for elt in size])

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