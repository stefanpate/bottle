from src import constants

from rdkit import Chem


def test_common_smarts():
    assert isinstance(constants.BaseSMARTS.CARBON, constants.BaseSMARTS)
    assert isinstance(constants.MoleculeSMARTS.C4_FG123.as_mol(), Chem.Mol)
