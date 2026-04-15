"""OPAM2-refined MolGpKa primitives.

A refinement of :mod:`src.molgpka` that incorporates upstream improvements
from the OPAM2 fork of MolGpKa. The original module is preserved unchanged;
this file provides new callables that can be used side-by-side.

Refinements over ``molgpka.py``:

* ``load_model`` is ``@lru_cache``-ed and uses ``weights_only=True``
  to avoid re-deserializing ``.pth`` files on every call.
* Prediction uses :class:`torch_geometric.nn.AttentionalAggregation` when
  available (the current PyG API), falling back to the deprecated
  :class:`GlobalAttention` for older installs.
* :func:`prepare_mol_for_prediction` replaces the SMILES round-trip in
  ``predict_for_protonate``: it preserves heavy-atom indices so callers
  can compare predictions against input atom indices.
* :func:`modify_acid` guards against ``GetNumExplicitHs() == 0`` so
  deprotonating atoms with no explicit hydrogens does not crash.
* :func:`cap_unstable_sites` caps combinatorial enumeration at
  :data:`MAX_UNSTABLE_SITES` (default 8 ⇒ 256 states max) by promoting
  overflow borderline sites to stable based on which side of ``ph``
  their pKa falls on.
* Optional ModelSeed-trained weight paths
  (``MODELSEED_{ACID,BASE}_MODEL_PATH``) for biochemistry-focused use.
"""

from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import torch
import torch.nn.functional as F
from torch.nn import BatchNorm1d, Linear

try:
    from torch_geometric.nn import AttentionalAggregation as _AttPool
except ImportError:  # older PyG
    from torch_geometric.nn import GlobalAttention as _AttPool

from torch_geometric.data import Data
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from .molgpka import (
    GCNConv,
    SMARTS_FILE,
    get_atom_features,
    get_bond_pair,
    get_ionization_aid,
    modify_base,
    n_features,
    hidden,
)

# ---------------------------------------------------------------------------
# Module-level paths
# ---------------------------------------------------------------------------

_MODULE_DIR = Path(__file__).parent
BASE_MODEL_PATH = _MODULE_DIR / "models" / "weight_base.pth"
ACID_MODEL_PATH = _MODULE_DIR / "models" / "weight_acid.pth"
MODELSEED_BASE_MODEL_PATH = _MODULE_DIR / "models" / "weight_base_modelseed.pth"
MODELSEED_ACID_MODEL_PATH = _MODULE_DIR / "models" / "weight_acid_modelseed.pth"

# Soft cap on combinatorial enumeration of borderline (near-pH) sites.
MAX_UNSTABLE_SITES = 8


# ---------------------------------------------------------------------------
# GCN architecture (uses AttentionalAggregation where available)
# ---------------------------------------------------------------------------


class GCNNetV2(torch.nn.Module):
    """5-layer GCN for per-atom pKa prediction with attentional pooling.

    Equivalent to :class:`src.molgpka.GCNNet` but uses the current PyG
    ``AttentionalAggregation`` layer instead of the deprecated
    ``GlobalAttention``.
    """

    def __init__(self) -> None:
        super().__init__()
        self.conv1 = GCNConv(n_features, 1024, cached=False)
        self.bn1 = BatchNorm1d(1024)
        self.conv2 = GCNConv(1024, 512, cached=False)
        self.bn2 = BatchNorm1d(512)
        self.conv3 = GCNConv(512, 256, cached=False)
        self.bn3 = BatchNorm1d(256)
        self.conv4 = GCNConv(256, 512, cached=False)
        self.bn4 = BatchNorm1d(512)
        self.conv5 = GCNConv(512, 1024, cached=False)
        self.bn5 = BatchNorm1d(1024)

        self.att = _AttPool(Linear(hidden, 1))
        self.fc2 = Linear(1024, 128)
        self.fc3 = Linear(128, 16)
        self.fc4 = Linear(16, 1)

    def reset_parameters(self) -> None:
        for conv in (self.conv1, self.conv2, self.conv3, self.conv4, self.conv5):
            conv.reset_parameters()
        self.att.reset_parameters()
        self.fc2.reset_parameters()
        self.fc3.reset_parameters()
        self.fc4.reset_parameters()

    def forward(self, data: Data) -> torch.Tensor:
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = F.relu(self.conv1(x, edge_index)); x = self.bn1(x)
        x = F.relu(self.conv2(x, edge_index)); x = self.bn2(x)
        x = F.relu(self.conv3(x, edge_index)); x = self.bn3(x)
        x = F.relu(self.conv4(x, edge_index)); x = self.bn4(x)
        x = F.relu(self.conv5(x, edge_index)); x = self.bn5(x)
        x = self.att(x, batch)
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        return self.fc4(x)


# ---------------------------------------------------------------------------
# Featurisation & model loading
# ---------------------------------------------------------------------------


def mol2vec(mol: Chem.Mol, atom_idx: int) -> Data:
    """Convert *mol* and ionisation site *atom_idx* into a PyG ``Data`` object.

    Evaluation-only (no pKa target).
    """
    import numpy as np

    node_f = get_atom_features(mol, atom_idx)
    edge_index = get_bond_pair(mol)
    batch = np.zeros(len(node_f))
    return Data(
        x=torch.tensor(node_f, dtype=torch.float32),
        edge_index=torch.tensor(edge_index, dtype=torch.long),
        batch=torch.tensor(batch, dtype=torch.long),
    )


@lru_cache(maxsize=4)
def load_model(model_file: str, device: str = "cpu") -> GCNNetV2:
    """Load a ``GCNNetV2`` checkpoint, cached by ``(path, device)``.

    Uses ``torch.load(..., weights_only=True)`` for safe deserialisation.
    """
    model = GCNNetV2().to(device)
    model.load_state_dict(
        torch.load(model_file, map_location=device, weights_only=True)
    )
    model.eval()
    return model


def model_pred(mol: Chem.Mol, aid: int, model: GCNNetV2, device: str = "cpu") -> float:
    """Run GCN inference and return the predicted pKa for atom *aid* in *mol*."""
    data = mol2vec(mol, aid).to(device)
    with torch.no_grad():
        pka_tensor = model(data)
    return float(pka_tensor.cpu().numpy()[0][0])


def predict_acid(
    mol: Chem.Mol,
    acid_model: str | Path | None = None,
    device: str = "cpu",
) -> dict[int, float]:
    """Return ``{atom_idx: predicted_pKa}`` for every acidic site in *mol*."""
    model = load_model(str(acid_model or ACID_MODEL_PATH), device)
    return {
        aid: model_pred(mol, aid, model, device)
        for aid in get_ionization_aid(mol, acid_or_base="acid")
    }


def predict_base(
    mol: Chem.Mol,
    base_model: str | Path | None = None,
    device: str = "cpu",
) -> dict[int, float]:
    """Return ``{atom_idx: predicted_pKa}`` for every basic site in *mol*."""
    model = load_model(str(base_model or BASE_MODEL_PATH), device)
    return {
        aid: model_pred(mol, aid, model, device)
        for aid in get_ionization_aid(mol, acid_or_base="base")
    }


def prepare_mol_for_prediction(mol: Chem.Mol, uncharged: bool = True) -> Chem.Mol:
    """Normalise *mol* for prediction while preserving heavy-atom indices.

    The older ``predict_for_protonate`` did a SMILES round-trip that
    canonicalised atom order and broke index-based comparisons.
    ``Uncharger`` alone preserves atom order, and ``AddHs`` appends
    hydrogens without renumbering heavy atoms, so heavy-atom indices in
    the returned mol match those of the input.
    """
    if uncharged:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    return AllChem.AddHs(mol)


def predict(
    mol: Chem.Mol,
    acid_model: str | Path | None = None,
    base_model: str | Path | None = None,
    uncharged: bool = True,
    device: str = "cpu",
) -> tuple[dict[int, float], dict[int, float]]:
    """Predict pKas for all acidic and basic sites in *mol*.

    Returns ``(base_dict, acid_dict)``.
    """
    mol = prepare_mol_for_prediction(mol, uncharged)
    return (
        predict_base(mol, base_model, device),
        predict_acid(mol, acid_model, device),
    )


def predict_for_protonate(
    mol: Chem.Mol,
    acid_model: str | Path | None = None,
    base_model: str | Path | None = None,
    uncharged: bool = True,
    device: str = "cpu",
) -> tuple[dict[int, float], dict[int, float], Chem.Mol]:
    """Like :func:`predict`, but also returns the H-added molecule used for
    prediction so callers can mutate it during protonation.
    """
    mol = prepare_mol_for_prediction(mol, uncharged)
    base_dict = predict_base(mol, base_model, device)
    acid_dict = predict_acid(mol, acid_model, device)
    return base_dict, acid_dict, mol


# ---------------------------------------------------------------------------
# Protonation-state manipulation
# ---------------------------------------------------------------------------


def modify_acid(at: Chem.Atom) -> None:
    """Set *at* to its deprotonated form, guarding against ``hnum == 0``."""
    hnum = at.GetNumExplicitHs()
    at.SetFormalCharge(-1)
    if hnum > 0:
        at.SetNumExplicitHs(hnum - 1)


def modify_mol(
    mol: Chem.Mol,
    acid_dict: dict[int, float],
    base_dict: dict[int, float],
) -> Chem.Mol:
    """Tag each atom in *mol* with its ionisation role (``"A"`` / ``"B"`` /
    ``"O"``) and predicted pKa, then strip explicit Hs."""
    for at in mol.GetAtoms():
        idx = at.GetIdx()
        if idx in acid_dict:
            # Acidic SMARTS centres point at the proton; charge its heavy neighbour.
            neighbor = at.GetNeighbors()[0]
            neighbor.SetProp("ionization", "A")
            neighbor.SetProp("pKa", str(acid_dict[idx]))
        elif idx in base_dict:
            at.SetProp("ionization", "B")
            at.SetProp("pKa", str(base_dict[idx]))
        else:
            at.SetProp("ionization", "O")
    return AllChem.RemoveHs(mol)


def get_pKa_data(
    mol: Chem.Mol,
    ph: float,
    tph: float,
) -> tuple[list, list, dict[int, float]]:
    """Partition tagged atoms into stable and unstable (borderline) sets."""
    stable_data: list = []
    unstable_data: list = []
    pkas: dict[int, float] = {}
    for at in mol.GetAtoms():
        props = at.GetPropsAsDict()
        acid_or_basic = props.get("ionization", False)
        if acid_or_basic not in ("A", "B"):
            continue
        pka = float(props.get("pKa", float("nan")))
        idx = at.GetIdx()
        pkas[idx] = pka
        if acid_or_basic == "A":
            if pka < ph - tph:
                stable_data.append([idx, pka, "A"])
            elif ph - tph <= pka <= ph + tph:
                unstable_data.append([idx, pka, "A"])
        else:
            if pka > ph + tph:
                stable_data.append([idx, pka, "B"])
            elif ph - tph <= pka <= ph + tph:
                unstable_data.append([idx, pka, "B"])
    return stable_data, unstable_data, pkas


def modify_stable_pka(new_mol: Chem.Mol, stable_data: list) -> None:
    """Apply protonation to all stable ionisation sites in *new_mol*."""
    for idx, _pka, acid_or_basic in stable_data:
        at = new_mol.GetAtomWithIdx(idx)
        if acid_or_basic == "A":
            modify_acid(at)
        elif acid_or_basic == "B":
            modify_base(at)


def modify_unstable_pka(
    mol: Chem.Mol,
    unstable_data: list,
    k: int,
) -> list[str]:
    """Enumerate all choices of *k* borderline sites to (de)protonate.

    Returns SMILES strings (one per combination).
    """
    new_smis: list[str] = []
    for pka_datas in combinations(unstable_data, k):
        if not pka_datas:
            continue
        new_mol = deepcopy(mol)
        for idx, _pka, acid_or_basic in pka_datas:
            at = new_mol.GetAtomWithIdx(idx)
            if acid_or_basic == "A":
                modify_acid(at)
            elif acid_or_basic == "B":
                modify_base(at)
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))
        new_smis.append(smi)
    return new_smis


def cap_unstable_sites(
    unstable_data: list,
    ph: float,
    max_sites: int = MAX_UNSTABLE_SITES,
) -> tuple[list, list]:
    """Cap combinatorial enumeration at *max_sites* borderline atoms.

    Overflow sites are promoted to stable based on which side of *ph*
    their pKa falls on; this keeps ``2**len(unstable)`` bounded for
    highly ionisable molecules.

    Returns ``(kept, promoted_stable)``.
    """
    if len(unstable_data) <= max_sites:
        return unstable_data, []
    sorted_sites = sorted(unstable_data, key=lambda d: abs(d[1] - ph))
    kept = sorted_sites[:max_sites]
    overflow = sorted_sites[max_sites:]
    promoted: list = []
    for idx, pka, acid_or_basic in overflow:
        if acid_or_basic == "A" and pka < ph:
            promoted.append([idx, pka, "A"])
        elif acid_or_basic == "B" and pka > ph:
            promoted.append([idx, pka, "B"])
    return kept, promoted


def _mol_to_format(mol: Chem.Mol, fmt: str):
    """Convert an RDKit ``Mol`` to the requested output format."""
    if fmt == "inchi":
        return Chem.MolToInchi(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    if fmt == "smiles":
        return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    if fmt == "mol":
        return mol
    raise ValueError(
        f"output_format must be one of 'inchi', 'smiles', 'mol'; got {fmt!r}"
    )


def opam_protonate_mol(
    smi_or_inchi: str,
    ph: float,
    tph: float,
    acid_model: str | Path | None = None,
    base_model: str | Path | None = None,
    output_format: str = "inchi",
    max_unstable_sites: int = MAX_UNSTABLE_SITES,
    device: str = "cpu",
) -> list:
    """Enumerate protonation microstates at ``ph ± tph``.

    Accepts a SMILES **or** an InChI string. Neutralises the input before
    pKa prediction so pre-charged zwitterions are handled correctly, then
    re-assigns charges based on predicted pKa vs target pH.

    Caps combinatorial expansion at *max_unstable_sites* borderline sites
    (default 8 → 256 states) via :func:`cap_unstable_sites`.
    """
    omol = Chem.MolFromSmiles(smi_or_inchi)
    if omol is None:
        omol = Chem.MolFromInchi(smi_or_inchi)
    if omol is None:
        raise ValueError(f"Could not parse input as SMILES or InChI: {smi_or_inchi!r}")

    base_dict, acid_dict, omol = predict_for_protonate(
        omol, acid_model=acid_model, base_model=base_model,
        uncharged=True, device=device,
    )
    mc = modify_mol(omol, acid_dict, base_dict)
    stable_data, unstable_data, _pkas = get_pKa_data(mc, ph, tph)

    unstable_data, promoted = cap_unstable_sites(unstable_data, ph, max_unstable_sites)
    stable_data = stable_data + promoted

    outputs: list = []
    n = len(unstable_data)
    if n == 0:
        new_mol = deepcopy(mc)
        modify_stable_pka(new_mol, stable_data)
        outputs.append(_mol_to_format(new_mol, output_format))
    else:
        for k in range(n + 1):
            new_mol = deepcopy(mc)
            modify_stable_pka(new_mol, stable_data)
            new_unsmis = modify_unstable_pka(new_mol, unstable_data, k)
            if output_format == "inchi":
                outputs.extend(Chem.MolToInchi(Chem.MolFromSmiles(s)) for s in new_unsmis)
            elif output_format == "smiles":
                outputs.extend(new_unsmis)
            elif output_format == "mol":
                outputs.extend(Chem.MolFromSmiles(s) for s in new_unsmis)
            else:
                raise ValueError(
                    f"output_format must be 'inchi', 'smiles', or 'mol'; got {output_format!r}"
                )
    return outputs
