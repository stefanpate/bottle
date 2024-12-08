from rdkit import Chem
from typing import Iterable
import pandas as pd
from src.cheminfo_utils import standardize_smiles
from collections import Counter
from src.operator_mapping import get_patts_from_operator_side
from src.utils import load_json

def generate_multisubstrate_rules(
        helper_molecules: Iterable[str],
        rules: pd.DataFrame,
        coreactants: pd.DataFrame,
        known_reactions: dict = {}
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Appends to rules and coreactants multisubstrate
    rules (and associated coreactants) for which one
    of helper molecules is the other substrate

    Args
    ----
    helper_molecules:Iterable[str]
        SMILES
    rules:pd.DataFrame
    coreactants:pd.DataFrame
    known_reactions:dict
        Known reactions mapped to rules with SMARTS aligned
        to rule template. If provided, cross-references multisubstrate
        rules against known reactions

    Returns
    -------
    rules:pd.DataFrame
    coreactants:pd.DataFrame
    '''

    if known_reactions:
        known_roles = _get_known_roles(known_reactions=known_reactions, rule_key='imt_rules')
    else:
        known_roles = {}

    helper_smiles = [standardize_smiles(smiles) for smiles in helper_molecules]
    helper_molecules = [Chem.MolFromSmiles(smiles) for smiles in helper_smiles]
    multi_any = rules[
        rules['Reactants'].apply(lambda x : Counter(x.split(';'))['Any'] > 1)
    ]
    tmp = []
    for hsmi, hmol in zip(helper_smiles, helper_molecules):
        for i, row in multi_any.iterrows():
            rule = row['Name']
            sma = row['SMARTS']
            patts = [Chem.MolFromSmarts(patt) for patt in get_patts_from_operator_side(sma, side=0)]
            for i, patt in enumerate(patts):
                if not hmol.HasSubstructMatch(patt):
                    continue
                    
                if not known_roles:
                    tmp.append(rule)
                    
                elif rule in known_roles and hsmi in known_roles[rule][i]:
                    tmp.append(rule)

    print()


def _get_known_roles(known_reactions: dict, rule_key: str) -> dict:
    '''
    Converts reaction dict with rule names under rule key
    and rule-aligned smarts under 'smarts'
    into dict indexed by rule name pointing to lists of 
    length = # reactants with set of molecules known to
    play role of operator reactant i at ith position
    '''
    tmp = {}
    for _, rxn in known_reactions.items():
        rcts = rxn['smarts'].split('>>')[0].split('.')
        # rcts = [standardize_smiles(smi) for smi in rxn['smarts'].split('>>')[0].split('.')] # NOTE Add back later
        if rxn[rule_key]:
            for rule in rxn[rule_key]:
                if rule not in tmp:
                    tmp[rule] = [[] for i in range(len(rcts))]

                for i, smi in enumerate(rcts):
                    tmp[rule][i].append(smi)

    known_roles = {}
    for k, v in tmp.items():
        known_roles[k] = tuple([set(elt) for elt in v])

    return known_roles




if __name__ == '__main__':
    from src.config import filepaths

    rules = pd.read_csv(
        filepath_or_buffer=filepaths['rules'] / 'JN3604IMT_rules.tsv',
        sep='\t'
    )

    coreactants = pd.read_csv(
        filepath_or_buffer=filepaths['coreactants'] / 'metacyc_coreactants.tsv',
        sep='\t'
    )

    starters = pd.read_csv(
        filepath_or_buffer=filepaths['starters_targets'] / 'ccm_v0.csv',
        sep=','
    )

    helper_molecules = list(starters['smiles'])

    succinate = "OC(=O)CCC(=O)O"
    pyruvate = 'CC(=O)C(=O)[O-]'

    known_reactions = load_json(filepaths['data'] / 'sprhea' / 'sprhea_240310_v3_mapped_no_subunits.json')


    generate_multisubstrate_rules(
        helper_molecules=helper_molecules,
        rules=rules,
        coreactants=coreactants,
        known_reactions=known_reactions
    )

    foo = {
        k: v 
        for k, v in known_reactions.items()
        if standardize_smiles(succinate) in v['smarts'].split('>>')[0].split('.')
        and standardize_smiles(pyruvate) in v['smarts'].split('>>')[0].split('.')  
        and v['imt_rules']
    }

print()