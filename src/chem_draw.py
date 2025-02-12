from rdkit import Chem
from rdkit.Chem import Draw
import svgutils.transform as st
from svgutils.compose import Unit
from collections import Counter
from src.config import filepaths
import numpy as np

plus = st.fromfile(filepaths['artifacts'] / 'plus.svg').to_str().decode()
arrow = st.fromfile(filepaths['artifacts'] / 'arrow.svg').to_str().decode()
double_border = st.fromfile(filepaths['artifacts'] / 'double_border.svg').to_str().decode()

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

# Temp

def _draw_reaction(pr_x_kr_img: st.SVGFigure, movex : int, rxn_sma: str, auto_scl: bool = False, size: tuple[int] = None) -> st.SVGFigure:
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

    if size:
        nrcts = len(reactants)
        npdts = len(products)
        arrow_width = int(st.fromstring(arrow).width.strip('px'))
        plus_width = int(st.fromstring(plus).width.strip('px'))
        tot_mol_width = size[0] - arrow_width - (nrcts + npdts - 2) * plus_width

        denominator = 0
        widths = ([], [])
        for i, side in enumerate([reactants, products]):
            for sub in side:
                mol = Chem.MolFromSmiles(sub)

                # Catch failed MolFromSmiles
                if mol is None: 
                    mol = Chem.MolFromSmiles(sub, sanitize=False)

                scl = max(np.log10(mol.GetNumAtoms()), 0.2)
                denominator += scl
                widths[i].append(scl)

        for j, side in enumerate(widths):
            for i, elt in enumerate(side):
                widths[j][i] = tot_mol_width * (elt / denominator)

    
    pr_x_kr_img, movex = _draw_side(reactants, movex, pr_x_kr_img, plus, auto_scl, widths[0])
    pr_x_kr_img, movex = _add_elt(arrow, movex, pr_x_kr_img)
    pr_x_kr_img, movex = _draw_side(products, movex, pr_x_kr_img, plus, auto_scl, widths[1])
    

    return pr_x_kr_img, movex

def _draw_pr_x_kr(pr: str, kr: str, size: tuple[int]):
    movex = 0
    pr_x_kr_img = st.SVGFigure()

    pr_x_kr_img, movex = _draw_reaction(pr_x_kr_img=pr_x_kr_img, movex=movex, rxn_sma=pr, size=(size[0] / 2, size[1]))
    pr_x_kr_img, movex = _add_elt(double_border, movex, pr_x_kr_img)
    pr_x_kr_img, movex = _draw_reaction(pr_x_kr_img=pr_x_kr_img, movex=movex, rxn_sma=kr, size=(size[0] / 2, size[1]))

    elt = st.fromstring(arrow)
    height = float(elt.height.strip('px')) # Assumed same for all elements
    pr_x_kr_img.width = Unit(movex)
    pr_x_kr_img.height = Unit(height)

    return pr_x_kr_img

if __name__ == '__main__':
    smi = 'CC(=O)O'
    mimg = draw_molecule(smi)

    rxn_sma = 'CCC.CCC.CCO>>CCCCCCCCO'
    fig = draw_reaction(rxn_sma, plus, arrow)

    _draw_pr_x_kr(pr=rxn_sma, kr=rxn_sma, size=(1500, 200))