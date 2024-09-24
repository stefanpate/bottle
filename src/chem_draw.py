from functools import partial
from typing import Callable
import PIL
from rdkit import Chem
from rdkit.Chem import Draw
import svgutils.transform as st
import svgutils.compose as sc
from collections import Counter
import numpy as np
import os

def draw_rxn_svg(rxn_sma, rid, hilite_atoms=None):
    fn = f"../artifacts/imgs/rxns/{rid}.svg"
    if os.path.exists(fn):
        pass
    else:
        reactants, products = [Counter(elt.split('.')) for elt in rxn_sma.split('>>')]

        elements = []
        movex = [0]
        element_ctr = 0
        elements, element_ctr, movex = draw_side(reactants, elements, element_ctr, movex, hilite_atoms=hilite_atoms)    
        elements.append(sc.SVG('../artifacts/arrow.svg'))
        movex.append(movex[element_ctr] + 40)
        element_ctr += 1
        elements, element_ctr, movex = draw_side(products, elements, element_ctr, movex, hilite_atoms=hilite_atoms)

        for i, elt in enumerate(elements):
            elt.moveto(movex[i], 0)

        rxn = sc.Figure(movex[-1], 200, *elements)
        rxn.save(fn)

    return fn


def draw_side(smi_stoich, elements, element_ctr, movex, hilite_atoms=None):
    mol_ctr = 0
    plus_ctr = 0
    for smi, stoich in smi_stoich.items():
        if hilite_atoms:
            fn = draw_mol_svg(smi, stoich, hilite_atoms[mol_ctr])
        else:
            fn = draw_mol_svg(smi, stoich)
        
        svg = sc.SVG(fn)
        width = svg.width
        elements.append(svg)
        movex.append(movex[element_ctr] + width)
        element_ctr +=1

        if plus_ctr < len(smi_stoich.values()) - 1:
            elements.append(sc.SVG('../artifacts/plus.svg'))
            movex.append(movex[element_ctr] + 40)
            element_ctr +=1

        plus_ctr += 1
        mol_ctr += 1

    return elements, element_ctr, movex

def draw_mol_svg(smiles, stoich, hilite_atoms=None):
    fn = f"../artifacts/imgs/mols/{hash((smiles, stoich))}.svg"
    if os.path.exists(fn):
        pass
    else:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None: # Catch MolFromSmiles troubles
            mol = Chem.MolFromSmiles(smiles, sanitize=False)

        nb = mol.GetNumAtoms()
        width = int(np.log10(nb) * 200) + 25
        d2d = Draw.MolDraw2DSVG(width, 200)

        if (stoich > 1) & (hilite_atoms is None):
            d2d.DrawMolecule(mol, legend=f"({stoich})", highlightAtoms=hilite_atoms)
        elif stoich > 1:
            d2d.DrawMolecule(mol, legend=f"({stoich})")
        elif (stoich == 1) & (hilite_atoms is None):
            d2d.DrawMolecule(mol)
        else:
            d2d.DrawMolecule(mol, highlightAtoms=hilite_atoms)
        
        d2d.FinishDrawing()
        fig = st.fromstring(d2d.GetDrawingText())
        fig.save(fn)
    
    return fn

# TODO PK added this code before discovering http://rdkit.org/docs/source/rdkit.Chem.PandasTools.html
# remove if this truly becomes dead code - otherwise, a decent pattern

DrawerProvider = Callable[[Chem.Mol], Draw.MolDraw2D]


def make_svg_drawer_non_scaling(width: int, height: int) -> DrawerProvider:
    def drawer_2d_non_scaling(_mol: Chem.Mol) -> Draw.MolDraw2D:
        return Draw.MolDraw2DSVG(width, height)
    return drawer_2d_non_scaling


def mol_to_svg(mol: Chem.Mol, drawer_provider: DrawerProvider) -> str:
    mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = drawer_provider(mol)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


df_formatter_mol_svg = partial(mol_to_svg, drawer_provider=make_svg_drawer_non_scaling(width=300, height=100))