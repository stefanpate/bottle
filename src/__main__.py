import sys
from typing import Annotated

import click
import typer
from rdkit import Chem

from src.smarts import MoleculeSmarts
from src.utils import filter_smiles_by_smarts, ensure_type

app = typer.Typer(no_args_is_help=True)

MolSmartsChoice = Annotated[
    str,
    typer.Option(click_type=click.Choice([mem.name for mem in MoleculeSmarts]))
]


@app.command()
def filter_smiles(
        in_smiles: Annotated[typer.FileText, typer.Argument()] = sys.stdin,
        out_smiles: Annotated[typer.FileTextWrite, typer.Argument()] = sys.stdout,
        *,
        smarts: str = None,
        mol_smarts: MolSmartsChoice = None,
):
    if not (smarts := smarts or (MoleculeSmarts[mol_smarts] if mol_smarts else None)):
        raise typer.BadParameter("You must provide exactly one of --smarts or --mol_smarts.")

    smarts = ensure_type(smarts, Chem.Mol, factory=Chem.MolFromSmarts)
    for line_num, line in enumerate(in_smiles, start=1):
        if not (smiles := ensure_type(line.strip(), Chem.Mol, factory=Chem.MolFromSmiles)):
            print(f'could not parse smiles on line {line_num}', file=sys.stderr)
            continue

        if smiles.HasSubstructMatch(smarts):
            out_smiles.write(line)


@app.callback()
def callback():
    """Dummy callback function. Remove when >1 subcommands."""
    pass


if __name__ == "__main__":
    app()
