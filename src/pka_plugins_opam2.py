"""OPAM2-refined ML-based pKa predictors for the eQuilibrator pipeline.

A refinement of :mod:`src.pka_plugins` that targets the OPAM2 fork of
MolGpKa. The original module is preserved unchanged; this file provides a
new :class:`MolGPKAOpam2` predictor that can be used side-by-side.

Key differences vs :class:`src.pka_plugins.MolGPKA`:

* Backed by :mod:`src.molgpka_opam2`, so it benefits from:
  - cached, ``weights_only=True`` model loading,
  - ``AttentionalAggregation`` pooling,
  - heavy-atom-index-preserving mol preparation,
  - safe deprotonation (``hnum == 0`` guard).
* Configurable acid/base model paths — defaults to the ModelSeed-trained
  weights when available, falling back to the MolGpKa stock weights.
* ``protonate_mol`` accepts a SMILES *or* InChI, a configurable output
  format (``"inchi"`` / ``"smiles"`` / ``"mol"``), and caps
  combinatorial enumeration via ``max_unstable_sites``.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from .pka_plugins import BasePkaPredictor
from .molgpka_opam2 import (
    ACID_MODEL_PATH,
    BASE_MODEL_PATH,
    MAX_UNSTABLE_SITES,
    MODELSEED_ACID_MODEL_PATH,
    MODELSEED_BASE_MODEL_PATH,
    cap_unstable_sites,
    get_pKa_data,
    load_model,
    model_pred,
    modify_mol,
    modify_stable_pka,
    modify_unstable_pka,
    opam_protonate_mol,
    predict_for_protonate,
)
from .molgpka import get_ionization_aid

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


def _default_model_path(preferred: Path, fallback: Path) -> Path:
    """Return *preferred* if it exists on disk, else *fallback*."""
    return preferred if Path(preferred).exists() else fallback


class MolGPKAOpam2(BasePkaPredictor):
    """OPAM2 pKa predictor: drop-in replacement for ChemAxon with OPAM2 fixes.

    Parameters
    ----------
    acid_model, base_model:
        Optional paths to custom acid/base weight files. Defaults to the
        ModelSeed-trained weights when available, otherwise the stock
        MolGpKa weights bundled under ``src/models/``.
    device:
        Torch device string (``"cpu"``, ``"cuda"``, ...).
    max_unstable_sites:
        Cap on combinatorial enumeration of borderline ionisation sites
        per molecule. ``2**max_unstable_sites`` protonation states at
        most.
    """

    def __init__(
        self,
        acid_model: str | Path | None = None,
        base_model: str | Path | None = None,
        device: str = "cpu",
        max_unstable_sites: int = MAX_UNSTABLE_SITES,
    ) -> None:
        self.acid_model_path = Path(
            acid_model
            or _default_model_path(MODELSEED_ACID_MODEL_PATH, ACID_MODEL_PATH)
        )
        self.base_model_path = Path(
            base_model
            or _default_model_path(MODELSEED_BASE_MODEL_PATH, BASE_MODEL_PATH)
        )
        self.device = device
        self.max_unstable_sites = max_unstable_sites
        # Warm the LRU cache so first prediction is not blocked by disk IO.
        self.acid_model = load_model(str(self.acid_model_path), device)
        self.base_model = load_model(str(self.base_model_path), device)

    # ------------------------------------------------------------------
    # BasePkaPredictor
    # ------------------------------------------------------------------

    def get_dissociation_constants(
        self,
        molecules: pd.DataFrame,
        error_log: str,
        num_acidic: int,
        num_basic: int,
        mid_ph: float,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Predict pKas and major microspecies for each row of *molecules*.

        See :meth:`src.pka_plugins.BasePkaPredictor.get_dissociation_constants`
        for full parameter documentation.
        """
        rows: list[dict] = []
        has_smiles_col = "smiles" in molecules.columns

        for row in molecules.itertuples():
            smi: str = (
                getattr(row, "smiles", None) if has_smiles_col else None
            ) or row.inchi

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                mol = Chem.MolFromInchi(smi) if smi else None
            if mol is None:
                logger.warning("Could not parse %r — skipping pKa prediction.", smi)
                rows.append({"major_ms": None})
                continue

            try:
                base_dict, acid_dict, _ = self.predict_for_protonate(mol)
            except Exception as exc:
                logger.warning("OPAM2 MolGPKA failed for %r: %s", smi, exc)
                rows.append({"major_ms": None})
                continue

            acid_pkas = sorted(acid_dict.values())[:num_acidic]
            base_pkas = sorted(base_dict.values(), reverse=True)[:num_basic]
            acid_pkas += [np.nan] * (num_acidic - len(acid_pkas))
            base_pkas += [np.nan] * (num_basic - len(base_pkas))

            row_dict: dict = {f"apKa{i + 1}": v for i, v in enumerate(acid_pkas)}
            row_dict.update({f"bpKa{i + 1}": v for i, v in enumerate(base_pkas)})
            row_dict["major_ms"] = self._get_major_microspecies(smi, mid_ph)
            rows.append(row_dict)

        results_df = pd.DataFrame(rows, index=molecules.index)
        view = molecules.join(results_df)

        pka_columns = (
            [f"apKa{i}" for i in range(1, num_acidic + 1)]
            + [f"bpKa{i}" for i in range(1, num_basic + 1)]
        )
        return view, pka_columns

    # ------------------------------------------------------------------
    # Prediction helpers
    # ------------------------------------------------------------------

    def predict_acid(self, mol: Chem.Mol) -> dict[int, float]:
        """Return ``{atom_idx: predicted_pKa}`` for every acidic site in *mol*."""
        return {
            aid: model_pred(mol, aid, self.acid_model, self.device)
            for aid in get_ionization_aid(mol, acid_or_base="acid")
        }

    def predict_base(self, mol: Chem.Mol) -> dict[int, float]:
        """Return ``{atom_idx: predicted_pKa}`` for every basic site in *mol*."""
        return {
            aid: model_pred(mol, aid, self.base_model, self.device)
            for aid in get_ionization_aid(mol, acid_or_base="base")
        }

    def predict_for_protonate(
        self,
        mol: Chem.Mol,
        uncharged: bool = True,
    ) -> tuple[dict[int, float], dict[int, float], Chem.Mol]:
        """Predict pKas for ionisable sites on the uncharged form of *mol*.

        Preserves heavy-atom indices (no SMILES round-trip). Returns
        ``(base_dict, acid_dict, mol_with_hs)``.
        """
        return predict_for_protonate(
            mol,
            acid_model=self.acid_model_path,
            base_model=self.base_model_path,
            uncharged=uncharged,
            device=self.device,
        )

    def protonate_mol(
        self,
        smi_or_inchi: str,
        ph: float,
        tph: float,
        output_format: str = "smiles",
        min_pka: float = 0.0,
        max_pka: float = 14.0,
    ) -> tuple[list, dict[int, float]]:
        """Enumerate protonation-state SMILES/InChIs/Mols for input at *ph* ± *tph*.

        Accepts SMILES or InChI. Caps combinatorial enumeration at
        ``self.max_unstable_sites`` borderline sites.

        Returns
        -------
        new_forms:
            List of enumerated protonation states in *output_format*.
        pkas:
            ``{atom_idx: pKa}`` for ionisable atoms within
            ``[min_pka, max_pka]``.
        """
        omol = Chem.MolFromSmiles(smi_or_inchi)
        if omol is None:
            omol = Chem.MolFromInchi(smi_or_inchi)
        if omol is None:
            logger.warning("Input %r could not be parsed by RDKit.", smi_or_inchi)
            return [], {}

        try:
            obase, oacid, omol_h = self.predict_for_protonate(omol)
        except Exception as exc:
            logger.warning("OPAM2 base/acid prediction failed for %r: %s", smi_or_inchi, exc)
            return [], {}

        mc = modify_mol(omol_h, oacid, obase)
        stable_data, unstable_data, pkas = get_pKa_data(mc, ph, tph)
        unstable_data, promoted = cap_unstable_sites(
            unstable_data, ph, self.max_unstable_sites
        )
        stable_data = stable_data + promoted

        outputs: list = []
        n = len(unstable_data)
        if n == 0:
            new_mol = deepcopy(mc)
            modify_stable_pka(new_mol, stable_data)
            smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))
            outputs.append(self._to_format(smi, output_format))
        else:
            for k in range(n + 1):
                new_mol = deepcopy(mc)
                modify_stable_pka(new_mol, stable_data)
                for smi in modify_unstable_pka(new_mol, unstable_data, k):
                    outputs.append(self._to_format(smi, output_format))

        pkas = {idx: pka for idx, pka in pkas.items() if min_pka < pka < max_pka}
        return outputs, pkas

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_format(smi: str, fmt: str):
        if fmt == "smiles":
            return smi
        if fmt == "inchi":
            return Chem.MolToInchi(Chem.MolFromSmiles(smi))
        if fmt == "mol":
            return Chem.MolFromSmiles(smi)
        raise ValueError(
            f"output_format must be 'smiles', 'inchi', or 'mol'; got {fmt!r}"
        )

    def _get_major_microspecies(self, smi_or_inchi: str, ph: float) -> str | None:
        """Return SMILES of the dominant protonation state at *ph*.

        Uses a ±0.5 pH tolerance window and returns the first enumerated
        form. Falls back to the input SMILES if enumeration produces no
        forms, or ``None`` if the input cannot be parsed.
        """
        mol = Chem.MolFromSmiles(smi_or_inchi) or Chem.MolFromInchi(smi_or_inchi)
        if mol is None:
            return None
        try:
            forms, _ = self.protonate_mol(
                smi_or_inchi, ph=ph, tph=0.5, output_format="smiles"
            )
            return forms[0] if forms else Chem.MolToSmiles(mol)
        except Exception as exc:
            logger.warning(
                "Major microspecies calculation failed for %r: %s", smi_or_inchi, exc
            )
            return None


# Convenience re-export so callers can do: from pka_plugins_opam2 import opam_protonate_mol
__all__ = ["MolGPKAOpam2", "opam_protonate_mol"]
