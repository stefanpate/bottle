from src.smarts import BaseSmarts, ActivationSmarts

from rdkit import Chem

import pytest


@pytest.mark.parametrize('smarts,smiles,has_substructure', [
    (ActivationSmarts.C4_FG123, 'CC(N)C(O)C(O)=O', True),
    (ActivationSmarts.C4_FG123, 'C(=O)(O)CC(O)C(O)=O', False),
    (ActivationSmarts.C4_FG123, 'OCCCCO', False),
    (ActivationSmarts.C4_FG123, 'CCCCCCCCC(N)C(O)C(O)=O', False),
    (ActivationSmarts.C4_FG123, 'CCCC', False),
])
def test_smiles_has_substruct(smarts: ActivationSmarts | str, smiles: str, has_substructure: bool):
    res = Chem.MolFromSmiles(smiles).HasSubstructMatch(smarts.as_mol())
    assert res == has_substructure


def test_sanity():
    assert isinstance(BaseSmarts.CARBON_ANY, BaseSmarts)
    assert isinstance(ActivationSmarts.C4_FG123.as_mol(), Chem.Mol)
