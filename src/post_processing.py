# Define classes for pathway and reaction entries
from src.utils import sort_x_by_y
from src.pathway_utils import get_stoich_pk
import PIL
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


class pathway:
    def __init__(self, rhashes, starter_hash=None, target_hash=None, prc_mcs=None, dG=None):
        self.starter = starter_hash
        self.target = target_hash
        self.rhashes = rhashes # Hash ids for the path's reactions, in order
        self.prc_mcs = prc_mcs # Peri-rxn-ctr MCS score ave over reactions
        self.dG = dG # Placeholder for thermo

class reaction:
    def __init__(self, rid, smarts, rules=[], known_rxns=[]):
        self.rid = rid
        self.smarts = smarts
        self.rules = rules
        self.known_rxns = known_rxns

    def sort_known_rxns(self):
        krs_w_mcs = [elt for elt in self.known_rxns if elt[0] is not None]
        krs_wo_mcs = [elt for elt in self.known_rxns if elt[0] is None]
        mcses = list(zip(*krs_w_mcs))[0]
        mean_mcses = list(map(lambda x: sum(x) / len(x), mcses))
        krs_w_mcs, _ = sort_x_by_y(krs_w_mcs, mean_mcses, reverse=True)
        self.known_rxns = list(krs_w_mcs) + krs_wo_mcs

def rxn_hash_2_rxn_sma(rhash, pk):
    '''
    Make reaction smarts string for
    reaction indexed by rhash in a pk
    object
    '''
    rxn_stoich = get_stoich_pk(rhash, pk)
    products = ".".join([".".join([smi]*stoich) for smi, stoich in rxn_stoich.items() if stoich >= 0])
    reactants = ".".join([".".join([smi]*abs(stoich)) for smi, stoich in rxn_stoich.items() if stoich <= 0])
    rxn_sma = ">>".join([reactants, products])
    return rxn_sma

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