from rdkit import Chem
from rdkit.Chem import Draw
import svgutils.transform as st
from svgutils.compose import Unit
from collections import Counter
from src.config import filepaths
import numpy as np

plus = st.fromfile(filepaths['artifacts'] / 'plus.svg')
arrow = st.fromfile(filepaths['artifacts'] / 'arrow.svg')

def draw_reaction(rxn_sma: str, plus: st.SVGFigure = plus, arrow: st.SVGFigure = arrow):
 
    reactants, products = [Counter(elt.split('.')) for elt in rxn_sma.split('>>')]

    movex = 0
    rxn_img = st.SVGFigure()
    def _add_elt(elt, movex, rxn_img):
        root = elt.getroot()
        root.moveto(movex, 0)
        rxn_img.append(root)
        movex += int(elt.width.strip('px'))
        return rxn_img, movex
    
    for i, (smi, stoich) in enumerate(reactants.items()):
        img = draw_molecule(smi, stoich)
        elt = st.fromstring(img)
        rxn_img, movex = _add_elt(elt, movex, rxn_img)

        if i < len(reactants) - 1:
            elt = plus
            rxn_img, movex = _add_elt(elt, movex, rxn_img)

    elt = arrow
    rxn_img, movex = _add_elt(elt, movex, rxn_img)

    for i, (smi, stoich) in enumerate(products.items()):
        img = draw_molecule(smi, stoich)
        elt = st.fromstring(img)
        rxn_img, movex = _add_elt(elt, movex, rxn_img)

        if i < len(products) - 1:
            elt = plus
            rxn_img, movex = _add_elt(elt, movex, rxn_img)

    height = float(elt.height.strip('px')) # Assumed same for all elements
    rxn_img.width = Unit(movex)
    rxn_img.height = Unit(height)

    rxn_img.save('test.svg')
    return rxn_img.to_str()

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