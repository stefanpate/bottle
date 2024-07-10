import sys
from enum import Enum

import typer

from src.smarts import MoleculeSMARTS
from src.utils import filter_smiles_by_smarts

app = typer.Typer()


class Boo(Enum):
    FOO = 'foozy'
    ZOO = 'zooby'


@app.command()
def filter_smiles(mol_smarts: Boo):
    print(mol_smarts)
    # pattern = Chem.MolFromSmarts(smarts)
    #
    # # Read SMILES from stdin
    # for line in sys.stdin:
    #     smiles = line.strip()
    #     mol = Chem.MolFromSmiles(smiles)
    #     if mol and mol.HasSubstructMatch(pattern):
    #         print(smiles)


if __name__ == "__main__":
    app()
