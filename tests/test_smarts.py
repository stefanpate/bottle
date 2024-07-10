from src.smarts import BaseSMARTS, MoleculeSMARTS

from rdkit import Chem

import pytest


@pytest.mark.parametrize('smarts,smiles,has_substructure', [
    (MoleculeSMARTS.C4_FG123, 'CC(N)C(O)C(O)=O', True),
    (MoleculeSMARTS.C4_FG123, 'C(=O)(O)CC(O)C(O)=O', False),
    (MoleculeSMARTS.C4_FG123, 'OCCCCO', False),
    (MoleculeSMARTS.C4_FG123, 'CCCCCCCCC(N)C(O)C(O)=O', False),
    (MoleculeSMARTS.C4_FG123, 'CCCC', False),
])
def test_smiles_has_substruct(smarts: MoleculeSMARTS | str, smiles: str, has_substructure: bool):
    res = Chem.MolFromSmiles(smiles).HasSubstructMatch(smarts.as_mol())
    assert res == has_substructure


def test_sanity():
    assert isinstance(BaseSMARTS.CARBON_ANY, BaseSMARTS)
    assert isinstance(MoleculeSMARTS.C4_FG123.as_mol(), Chem.Mol)
