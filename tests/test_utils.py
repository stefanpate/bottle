import pytest

from src.utils import filter_smiles_by_smarts
from src.smarts import ActivationSmarts


@pytest.mark.parametrize('smarts,smiles_in,smiles_out', [
    (
            ActivationSmarts.C4_FG123,
            ['CC(N)C(O)C(O)=O', 'C(=O)(O)CC(O)C(O)=O', 'OCCCCO', 'CCCCCCCCC(N)C(O)C(O)=O', 'CCCC'],
            ['CC(N)C(O)C(O)=O'],
    )
])
def test_filter_smiles_by_smarts(smarts: str, smiles_in: list[str], smiles_out: list[str]):
    res = list(filter_smiles_by_smarts(smiles_in, smarts))
    assert res == smiles_out
