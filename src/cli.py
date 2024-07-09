import sys
import typer
from src.constants import MoleculeSMARTS
from src.utils import filter_smiles_by_smarts

app = typer.Typer()


@app.command()
def filter_smiles(mol_smarts: MoleculeSMARTS | str):
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
