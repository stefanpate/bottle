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

def draw_rxn_svg(rxn_sma):
    reactants, products = [elt.split('.') for elt in rxn_sma.split('>>')]
    reactants, products = Counter(reactants), Counter(products)

    fns = []
    movex = [0]
    plus_ctr = 0
    element_ctr = 0
    for smi, stoich in reactants.items():
        fn, width = draw_mol_svg(smi, stoich)
        fns.append(fn)
        movex.append(movex[element_ctr] + width)
        element_ctr +=1

        if plus_ctr < len(reactants.values()) - 1:
            fns.append('../artifacts/plus.svg')
            movex.append(movex[element_ctr] + 40)
            element_ctr +=1

        plus_ctr += 1
        
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
    return rxn

def draw_mol_svg(smiles, stoich):
    mol = Chem.MolFromSmiles(smiles)
    nb = mol.GetNumAtoms()
    width = int(np.log10(nb) * 200) + 25
    d2d = Draw.MolDraw2DSVG(width, 200)

    if stoich > 1:
        d2d.DrawMolecule(mol, legend=f"({stoich})")
    else:
        d2d.DrawMolecule(mol)

    d2d.FinishDrawing()
    fn = f"../artifacts/mol_svgs/{hash((smiles, stoich))}.svg"
    fig = st.fromstring(d2d.GetDrawingText())
    fig.save(fn)
    return fn, width

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