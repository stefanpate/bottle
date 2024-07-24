import sys
from contextlib import nullcontext
from typing import Annotated

from cyclopts import App, Group, Parameter, validators
from cyclopts.types import File
from rdkit import Chem

from src.smarts import ActivationSmarts
from src.utils import ensure_type

app = App("Bottle CLI")

_smarts_group = Group(
    'smarts (choose one)',
    default_parameter=Parameter(negative=""),  # Disable "--no-" flags
    validator=validators.LimitedChoice(min=1),  # Mutually Exclusive Options
)


@app.command()
def filter_smiles(
        in_smiles: Annotated[File, Parameter(name=('--input', '-I'))] = None,
        *,
        out_smiles: Annotated[File, Parameter(name=('--output', '-O'))] = None,
        smarts: Annotated[str, Parameter(group=_smarts_group)] = None,
        mol_smarts: Annotated[ActivationSmarts, Parameter(group=_smarts_group)] = None,
):
    """Filter a stream of SMILES lines by a SMARTS pattern.

    Args:
        in_smiles: input stream of SMILES lines
        out_smiles: output stream of SMILES lines
        smarts: free-form SMARTS pattern
        mol_smarts: predefined SMARTS patterns
    """
    smarts_str = mol_smarts or smarts
    smarts = ensure_type(smarts_str, Chem.Mol, factory=Chem.MolFromSmarts)

    with (
        in_smiles.open(mode='r') if in_smiles else nullcontext(sys.stdin) as fin,
        out_smiles.open(mode='w') if out_smiles else nullcontext(sys.stdout) as fout,
    ):
        for line_num, line in enumerate(fin, start=1):
            if not (smiles := ensure_type(line.strip(), Chem.Mol, factory=Chem.MolFromSmiles)):
                print(f'could not parse smiles on line {line_num}', file=sys.stderr)
                continue

            if smiles.HasSubstructMatch(smarts):
                fout.write(line)


def main():
    app()


if __name__ == "__main__":
    main()
