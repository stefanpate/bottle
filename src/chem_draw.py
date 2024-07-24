from functools import partial
from typing import Callable

import PIL
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import svgutils.transform as st
import svgutils.compose as sc
from collections import Counter
import numpy as np

'''
SVG
'''

def draw_pwy_svg(sma_hash_pairs, path_id, pwy_fn=None):
    fns = []
    widths = []
    for pair in sma_hash_pairs:
        pred_fn, pred_width = draw_rxn_svg(*pair[0])
        known_fn, known_width = draw_rxn_svg(*pair[1])
        fns.append((pred_fn, known_fn))
        widths.append((pred_width, known_width))

    max_pred_width = max(list(zip(*widths))[0])
    max_known_width = max(list(zip(*widths))[1])
    elements = []
    for i, row in enumerate(fns):
        for j, half in enumerate(row):
            pred_delta = max_pred_width - widths[i][j]
            elements.append(sc.SVG(half).move(pred_delta + j * (max_pred_width + 40 - pred_delta), 200 * i))

        elements.append(sc.SVG('../artifacts/double_border.svg').move(max_pred_width, 200 * i))
    
    elements.append(sc.Text(f"{path_id}", 25, 25, size=16, weight='bold')) # Path id
    pwy_svg = sc.Figure(max_pred_width + 40 + max_known_width, 200 * len(fns),
            *elements
            )
    
    if pwy_fn:
        pwy_svg.save(pwy_fn)
    else:
        return pwy_svg
        
def draw_rxn_svg(rxn_sma, rhash=None, hilite_atoms=None):
    reactants, products = [elt.split('.') for elt in rxn_sma.split('>>')]
    reactants, products = Counter(reactants), Counter(products)

    fns = []
    movex = [0]
    reactant_ctr = 0
    plus_ctr = 0
    element_ctr = 0
    for smi, stoich in reactants.items():
        if hilite_atoms:
            fn, width = draw_mol_svg(smi, stoich, hilite_atoms[reactant_ctr])
        else:
            fn, width = draw_mol_svg(smi, stoich)
        fns.append(fn)
        movex.append(movex[element_ctr] + width)
        element_ctr +=1

        if plus_ctr < len(reactants.values()) - 1:
            fns.append('../artifacts/plus.svg')
            movex.append(movex[element_ctr] + 40)
            element_ctr +=1

        plus_ctr += 1
        reactant_ctr += 1
        
    fns.append('../artifacts/arrow.svg')
    movex.append(movex[element_ctr] + 40)
    element_ctr +=1


    plus_ctr = 0
    for smi, stoich in products.items():
        fn, width = draw_mol_svg(smi, stoich)
        fns.append(fn)
        movex.append(movex[element_ctr] + width)
        element_ctr +=1

        if plus_ctr < len(products.values()) - 1:
            fns.append('../artifacts/plus.svg')
            movex.append(movex[element_ctr] + 40)
            element_ctr +=1

        plus_ctr += 1

    elements = [sc.SVG(elt) for elt in fns]
    for i, elt in enumerate(elements):
        elt.moveto(movex[i], 0)

    rxn = sc.Figure(movex[-1], 200,
            *elements
            )
    width = movex[-1]

    if rhash:
        fn = f"../artifacts/rxn_svgs/{rhash}.svg"
        rxn.save(fn)

    return fn, width # TODO FIX. If not supplying rhash, will give strange result for fn

def draw_mol_svg(smiles, stoich, hilite_atoms=None):
    mol = Chem.MolFromSmiles(smiles)
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
    fn = f"../artifacts/mol_svgs/{hash((smiles, stoich))}.svg"
    fig = st.fromstring(d2d.GetDrawingText())
    fig.save(fn)
    return fn, width


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
    return drawer.GetDrawingText()


df_formatter_mol_svg = partial(mol_to_svg, drawer_provider=make_svg_drawer_non_scaling(width=300, height=100))


'''
PNG
'''

# Pathway drawing functions
def draw_rxn(rxn_sma):
    return Draw.ReactionToImage(
        AllChem.ReactionFromSmarts(rxn_sma, useSmiles=True),
        subImgSize=(200, 200), useSVG=False, drawOptions=None, returnPNG=False
    )

def get_concat_h(im1, im2):
    dst = PIL.Image.new('RGB', (im1.width + im2.width, max(im1.height, im2.height)))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    dst = PIL.Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def draw_pathway(pred_known_pairs):
    for i, elt in enumerate(pred_known_pairs):
        if i == 0:
            img = get_concat_h(*elt)
        else:
            img = get_concat_v(img, get_concat_h(*elt))

    return img