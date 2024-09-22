import collections
import json
import orjson
import pandas as pd
from pathlib import Path
from typing import Iterable, TypeVar

from rdkit import Chem
from rdkit.Chem import BRICS

T = TypeVar('T')


def save_json(data, save_to):
    with open(save_to, 'w') as f:
        json.dump(data, f)


def load_json(path):
    json_str = Path(path).read_text('utf-8')
    data = orjson.loads(json_str)
    return data


# def ensure_type(obj: Any, target_type: Type[T], *, factory: Callable[[Any], T] | None = None) -> T:
#     factory = factory or target_type
#     return obj if isinstance(obj, target_type) else factory(obj)


# def filter_smiles_by_smarts(
#         smiles: Iterable[str | rdkit.Chem.Mol],
#         smarts: str | rdkit.Chem.Mol
# ) -> Iterable[rdkit.Chem.Mol]:
#     smarts = ensure_type(smarts, rdkit.Chem.Mol, factory=rdkit.Chem.MolFromSmarts)
#     return (
#         elem for elem in smiles
#         if ensure_type(elem, rdkit.Chem.Mol, factory=rdkit.Chem.MolFromSmiles).HasSubstructMatch(smarts)
#     )


def BRICSDecompositionsToFrame(mols: Iterable[Chem.Mol], *, keep_decomp_tuple=True, **kwargs) -> pd.DataFrame:
    """BRICS Decomposition results to DataFrame. Especially useful when rendered with rdkit.Chem.PandasTools

    A convenience function for creating a DataFrame from a set of molecule BRICS decompositions.
    See http://www.rdkit.org/new_docs/source/rdkit.Chem.BRICS.html#module-rdkit.Chem.BRICS

    Args:
        mols: iterable of rdkit.Chem.Mol
        keep_decomp_tuple: preserve the column with the decomp result tuple
        **kwargs: keyword args to pass to BRICS.BRICSDecompose

    Returns:
        DataFrame with BRICS decompositions
    """
    sr_decomp = pd.Series(mols).apply(lambda mol: BRICS.BRICSDecompose(mol, returnMols=True, **kwargs))
    df = sr_decomp.apply(pd.Series).rename(columns=lambda idx: f'brics_frag_{idx + 1}')
    if keep_decomp_tuple:
        df.insert(0, 'brics_decomp', sr_decomp)
    return df
