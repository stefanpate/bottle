"""ML-based pKa predictors that can replace ChemAxon in the eQuilibrator pipeline.

Typical usage (monkey-patch before any calls to ``lc.add_compounds``):

    import equilibrator_assets.chemaxon as _chemaxon
    from src.pka_plugins import MolGPKA

    _predictor = MolGPKA()
    _chemaxon.get_dissociation_constants = _predictor.get_dissociation_constants
"""

import logging
from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path

import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from .molgpka import (
    BASE_MODEL_PATH,
    ACID_MODEL_PATH,
    MAX_UNSTABLE_SITES,
    MODELSEED_ACID_MODEL_PATH,
    MODELSEED_BASE_MODEL_PATH,
    cap_unstable_sites,
    get_ionization_aid,
    load_model,
    model_pred,
    modify_mol,
    modify_stable_pka,
    modify_unstable_pka,
    get_pKa_data,
    prepare_mol_for_prediction,
)

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


class BasePkaPredictor(ABC):
    """Abstract base class for pKa predictors compatible with the eQuilibrator pipeline."""

    @abstractmethod
    def get_dissociation_constants(
        self,
        molecules: pd.DataFrame,
        error_log: str,
        num_acidic: int,
        num_basic: int,
        mid_ph: float,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Compute dissociation constants and the major microspecies at *mid_ph*.

        Drop-in replacement for ``equilibrator_assets.chemaxon.get_dissociation_constants``.

        Parameters
        ----------
        molecules:
            DataFrame with at least an ``inchi`` column (may contain SMILES) and
            optionally a ``smiles`` column.  Matches the format produced by
            ``equilibrator_assets.generate_compound.create_compounds``.
        error_log:
            Base path for error output files (accepted for API compatibility; ML
            predictors may log via Python's logging instead).
        num_acidic:
            Maximum number of acidic pKa columns to return (``apKa1`` … ``apKaN``).
        num_basic:
            Maximum number of basic pKa columns to return (``bpKa1`` … ``bpKaN``).
        mid_ph:
            pH used to determine the major microspecies.

        Returns
        -------
        view : pd.DataFrame
            Input *molecules* DataFrame joined with pKa columns and ``major_ms``.
        pka_columns : list[str]
            Names of the pKa columns present in *view*
            (``["apKa1", ..., "bpKa1", ...]``).
        """
        ...


class MolGPKA(BasePkaPredictor):
    """pKa predictor backed by a graph-convolutional network (MolGPKA).

    Loads separate GCNNet weights for acidic and basic site prediction from
    ``src/models/``.
    """

    def __init__(self) -> None:
        self.base_model = load_model(BASE_MODEL_PATH)
        self.acid_model = load_model(ACID_MODEL_PATH)

    # ------------------------------------------------------------------
    # Public API (BasePkaPredictor)
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

        See ``BasePkaPredictor.get_dissociation_constants`` for full parameter
        documentation.
        """
        rows: list[dict] = []
        has_smiles_col = "smiles" in molecules.columns

        for row in molecules.itertuples():
            smi: str = (
                getattr(row, "smiles", None) if has_smiles_col else None
            ) or row.inchi

            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning("Could not parse SMILES %r — skipping pKa prediction.", smi)
                rows.append({"major_ms": None})
                continue

            try:
                base_dict, acid_dict, _ = self.predict_for_protonate(mol)
            except Exception as exc:
                logger.warning(
                    "MolGPKA prediction failed for %r: %s", smi, exc
                )
                rows.append({"major_ms": None})
                continue

            # Collect pKas: acids ascending (most acidic first), bases descending
            acid_pkas = sorted(acid_dict.values())[:num_acidic]
            base_pkas = sorted(base_dict.values(), reverse=True)[:num_basic]

            # Pad with nulls if there are fewer than num_acidic/basic sites
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
        """Return ``{atom_idx: predicted_pKa}`` for all acidic sites in *mol*."""
        return {
            aid: model_pred(mol, aid, self.acid_model)
            for aid in get_ionization_aid(mol, acid_or_base="acid")
        }

    def predict_base(self, mol: Chem.Mol) -> dict[int, float]:
        """Return ``{atom_idx: predicted_pKa}`` for all basic sites in *mol*."""
        return {
            aid: model_pred(mol, aid, self.base_model)
            for aid in get_ionization_aid(mol, acid_or_base="base")
        }

    def predict_for_protonate(
        self,
        mol: Chem.Mol,
        uncharged: bool = True,
    ) -> tuple[dict[int, float], dict[int, float], Chem.Mol]:
        """Predict pKas for all ionisable sites on the uncharged form of *mol*.

        Parameters
        ----------
        mol:
            Input molecule.
        uncharged:
            If True, neutralise formal charges before prediction.

        Returns
        -------
        base_dict, acid_dict, mol_with_h:
            Per-atom pKa dicts and the Hs-added molecule used for prediction.
        """
        if uncharged:
            un = rdMolStandardize.Uncharger()
            mol = un.uncharge(mol)
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        mol = AllChem.AddHs(mol)
        base_dict = self.predict_base(mol)
        acid_dict = self.predict_acid(mol)
        return base_dict, acid_dict, mol

    def protonate_mol(
        self,
        smi: str,
        ph: float,
        tph: float,
        min_pka: float = 0.0,
        max_pka: float = 14.0,
    ) -> tuple[list[str], dict[int, float]]:
        """Enumerate protonation-state SMILES for *smi* at *ph* ± *tph*.

        Parameters
        ----------
        smi:
            Input SMILES string.
        ph:
            Target pH.
        tph:
            Half-width of the "uncertain" pH window.  Atoms whose pKa falls
            within ``ph ± tph`` are treated as potentially in either state.
        min_pka, max_pka:
            pKa range filter applied to the returned *pkas* dict.

        Returns
        -------
        new_smis:
            List of SMILES for each enumerated protonation state.
        pkas:
            ``{atom_idx: pKa}`` for ionisable atoms within [min_pka, max_pka].
        """
        omol = Chem.MolFromSmiles(smi)
        if omol is None:
            logger.warning("SMILES %r could not be parsed by RDKit.", smi)
            return [], {}

        try:
            obase_dict, oacid_dict, omol = self.predict_for_protonate(omol)
        except Exception as exc:
            logger.warning(
                "MolGPKA base/acid prediction failed for %r: %s", smi, exc
            )
            return [], {}

        mc = modify_mol(omol, oacid_dict, obase_dict)
        stable_data, unstable_data, pkas = get_pKa_data(mc, ph, tph)

        new_smis: list[str] = []
        n = len(unstable_data)
        if n == 0:
            new_mol = deepcopy(mc)
            modify_stable_pka(new_mol, stable_data)
            new_smis.append(Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol))))
        else:
            for i in range(n + 1):
                new_mol = deepcopy(mc)
                modify_stable_pka(new_mol, stable_data)
                new_smis.extend(modify_unstable_pka(new_mol, unstable_data, i))

        pkas = {idx: pka for idx, pka in pkas.items() if min_pka < pka < max_pka}
        return new_smis, pkas

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_major_microspecies(self, smi: str, ph: float) -> str | None:
        """Return the SMILES of the dominant protonation state at *ph*.

        Uses ``protonate_mol`` with a ±0.5 pH tolerance window and returns the
        first enumerated form.  Falls back to the input SMILES if enumeration
        produces no forms, or ``None`` if the input cannot be parsed.
        """
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        try:
            major_smis, _ = self.protonate_mol(smi, ph=ph, tph=0.5)
            return major_smis[0] if major_smis else smi
        except Exception as exc:
            logger.warning(
                "Major microspecies calculation failed for %r: %s", smi, exc
            )
            return None


class MolGPKAOpam2(MolGPKA):
    """OPAM2-refined pKa predictor.

    Inherits from :class:`MolGPKA` — transitively from
    :class:`BasePkaPredictor` — so the eQuilibrator pipeline's drop-in
    contract is preserved and the stock MolGpKa flow (plain
    ``MolGPKA()``) continues to work unchanged.

    Differences vs :class:`MolGPKA`:

    * Defaults to the ModelSeed-trained weights when they are present on
      disk (``weight_{acid,base}_modelseed.pth``); falls back to the
      stock MolGpKa weights otherwise. Both paths may be overridden via
      the constructor.
    * :meth:`predict_for_protonate` uses
      :func:`src.molgpka.prepare_mol_for_prediction`, which preserves
      heavy-atom indices between the input and the H-added mol (no
      SMILES round-trip).
    * :meth:`protonate_mol` accepts a SMILES **or** an InChI, caps
      combinatorial enumeration via
      :func:`src.molgpka.cap_unstable_sites`, and supports
      ``output_format`` in ``{"smiles", "inchi", "mol"}``.
    * :meth:`get_dissociation_constants` falls back to
      ``Chem.MolFromInchi`` when SMILES parsing fails, so InChI-only
      rows are handled.
    """

    def __init__(
        self,
        acid_model: Path | str | None = None,
        base_model: Path | str | None = None,
        device: str = "cpu",
        max_unstable_sites: int = MAX_UNSTABLE_SITES,
    ) -> None:
        self.acid_model_path = Path(
            acid_model
            or (
                MODELSEED_ACID_MODEL_PATH
                if Path(MODELSEED_ACID_MODEL_PATH).exists()
                else ACID_MODEL_PATH
            )
        )
        self.base_model_path = Path(
            base_model
            or (
                MODELSEED_BASE_MODEL_PATH
                if Path(MODELSEED_BASE_MODEL_PATH).exists()
                else BASE_MODEL_PATH
            )
        )
        self.device = device
        self.max_unstable_sites = max_unstable_sites
        self.base_model = load_model(self.base_model_path, device)
        self.acid_model = load_model(self.acid_model_path, device)

    # ------------------------------------------------------------------
    # Overrides
    # ------------------------------------------------------------------

    def predict_acid(self, mol: Chem.Mol) -> dict[int, float]:
        return {
            aid: model_pred(mol, aid, self.acid_model, device=self.device)
            for aid in get_ionization_aid(mol, acid_or_base="acid")
        }

    def predict_base(self, mol: Chem.Mol) -> dict[int, float]:
        return {
            aid: model_pred(mol, aid, self.base_model, device=self.device)
            for aid in get_ionization_aid(mol, acid_or_base="base")
        }

    def predict_for_protonate(
        self,
        mol: Chem.Mol,
        uncharged: bool = True,
    ) -> tuple[dict[int, float], dict[int, float], Chem.Mol]:
        """Predict pKas while preserving heavy-atom indices of *mol*."""
        mol = prepare_mol_for_prediction(mol, uncharged=uncharged)
        return self.predict_base(mol), self.predict_acid(mol), mol

    def protonate_mol(
        self,
        smi: str,
        ph: float,
        tph: float,
        min_pka: float = 0.0,
        max_pka: float = 14.0,
        output_format: str = "smiles",
    ) -> tuple[list, dict[int, float]]:
        """Enumerate protonation states at ``ph ± tph``.

        Accepts SMILES or InChI. Caps combinatorial enumeration at
        ``self.max_unstable_sites``. Returns enumerated forms in
        *output_format* and the ``{atom_idx: pKa}`` map for ionisable
        atoms within ``[min_pka, max_pka]``.
        """
        omol = Chem.MolFromSmiles(smi)
        if omol is None:
            omol = Chem.MolFromInchi(smi)
        if omol is None:
            logger.warning("Input %r could not be parsed by RDKit.", smi)
            return [], {}

        try:
            obase_dict, oacid_dict, omol = self.predict_for_protonate(omol)
        except Exception as exc:
            logger.warning(
                "MolGPKAOpam2 base/acid prediction failed for %r: %s", smi, exc
            )
            return [], {}

        mc = modify_mol(omol, oacid_dict, obase_dict)
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
            smi_out = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)))
            outputs.append(self._to_format(smi_out, output_format))
        else:
            for k in range(n + 1):
                new_mol = deepcopy(mc)
                modify_stable_pka(new_mol, stable_data)
                for new_smi in modify_unstable_pka(new_mol, unstable_data, k):
                    outputs.append(self._to_format(new_smi, output_format))

        pkas = {idx: pka for idx, pka in pkas.items() if min_pka < pka < max_pka}
        return outputs, pkas

    def get_dissociation_constants(
        self,
        molecules: pd.DataFrame,
        error_log: str,
        num_acidic: int,
        num_basic: int,
        mid_ph: float,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Like :meth:`MolGPKA.get_dissociation_constants`, but also tolerates
        rows whose only available representation is an InChI."""
        rows: list[dict] = []
        has_smiles_col = "smiles" in molecules.columns
        for row in molecules.itertuples():
            smi: str = (
                getattr(row, "smiles", None) if has_smiles_col else None
            ) or row.inchi

            mol = Chem.MolFromSmiles(smi) if smi else None
            if mol is None and smi:
                mol = Chem.MolFromInchi(smi)
            if mol is None:
                logger.warning("Could not parse %r — skipping pKa prediction.", smi)
                rows.append({"major_ms": None})
                continue

            try:
                base_dict, acid_dict, _ = self.predict_for_protonate(mol)
            except Exception as exc:
                logger.warning("MolGPKAOpam2 prediction failed for %r: %s", smi, exc)
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
