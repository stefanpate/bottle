"""eQuilibrator integration helpers ported from OPAM2.

Thin wrapper around :mod:`equilibrator_api` for computing transformed
Gibbs free energies of formation (ΔG'°) of compounds after protonation
by the OPAM2 pipeline. The heavy lifting is done by ``equilibrator-api``,
which must be installed separately (it is an optional dependency).

Notes (eQuilibrator 0.6 API):

* No ``cc.get_transformed_dg(smiles)`` method exists. Instead:
    1. Convert SMILES → InChI → InChI-Key via RDKit.
    2. Neutralise the molecule first (``Uncharger``) so the 14-char
       InChI-Key skeleton prefix is protonation-invariant.
    3. Look up the compound via
       ``cc.search_compound_by_inchi_key(skeleton)``.
    4. Compute ``cc.standard_dg_prime(Reaction({compound: 1}))``.
* eQuilibrator canonicalises protonation internally and applies its own
  Legendre transform, so the protonation state of the input SMILES does
  not affect the ΔG'° output. The skeleton lookup absorbs any
  protonation difference.
"""

from __future__ import annotations

from typing import Mapping

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

_cc = None
_Q_ = None
_Reaction = None
_uncharger = rdMolStandardize.Uncharger()


def _ensure_equilibrator():
    """Lazily load the ``ComponentContribution`` model and unit registry.

    Caches the loaded objects in module globals so subsequent calls are
    free.
    """
    global _cc, _Q_, _Reaction
    if _cc is not None:
        return _cc, _Q_, _Reaction
    try:
        from equilibrator_api import ComponentContribution, Q_, Reaction
    except ImportError as exc:
        raise ImportError(
            "equilibrator-api is required for calculate_delta_g but is not "
            "installed. Install it with `pip install equilibrator-api`."
        ) from exc
    _cc = ComponentContribution()
    _Q_ = Q_
    _Reaction = Reaction
    return _cc, _Q_, _Reaction


def _smiles_to_skeleton_key(smiles: str) -> str | None:
    """Convert a SMILES to the 14-char InChI-Key skeleton prefix.

    Neutralises formal charges first so the skeleton is
    protonation-invariant.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = _uncharger.uncharge(mol)
    if mol is None:
        return None
    inchi = Chem.MolToInchi(mol)
    if not inchi:
        return None
    inchi_key = Chem.InchiToInchiKey(inchi)
    if not inchi_key or len(inchi_key) < 14:
        return None
    return inchi_key[:14]


def _resolve_compound(cc, smiles: str, compound_id: str):
    """Resolve a SMILES to an eQuilibrator Compound via skeleton key."""
    skeleton = _smiles_to_skeleton_key(smiles)
    if skeleton is None:
        return None
    hits = cc.search_compound_by_inchi_key(skeleton)
    if not hits:
        return None
    return hits[0]


def _coerce_to_smiles(
    structures: Mapping[str, str | None],
    compound_id: str,
) -> str:
    """Pick the best available structure representation and return SMILES.

    Priority: SMILES > MOL > InChI.
    """
    smi = structures.get("smiles")
    if smi:
        if Chem.MolFromSmiles(smi) is None:
            raise ValueError(f"Invalid SMILES for {compound_id}: {smi!r}")
        return smi

    mol_block = structures.get("mol")
    if mol_block:
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            raise ValueError(f"Invalid MOL block for {compound_id}")
        return Chem.MolToSmiles(mol)

    inchi = structures.get("inchi")
    if inchi:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            raise ValueError(f"Invalid InChI for {compound_id}: {inchi!r}")
        return Chem.MolToSmiles(mol)

    raise ValueError(f"No structure representation provided for {compound_id}")


def calculate_delta_g(
    hash_of_structures: Mapping[str, Mapping[str, str | None]],
    pH: float = 7.0,
    ionic_strength_M: float = 0.25,
) -> dict[str, float | None]:
    """Compute transformed ΔG'° (kJ/mol) for each compound.

    Parameters
    ----------
    hash_of_structures
        Mapping from compound ID to a dict of structure representations.
        Each inner dict may contain any subset of ``"smiles"``,
        ``"mol"``, ``"inchi"``. The first non-empty representation is
        used, in that priority order.
    pH
        Target pH for the transformed ΔG calculation. Default 7.0.
    ionic_strength_M
        Ionic strength in molar. Default 0.25.

    Returns
    -------
    dict[str, float | None]
        Compound ID → ΔG'° in kJ/mol, or ``None`` if the calculation
        failed for that compound (with the error printed to stdout).
    """
    cc, Q_, Reaction = _ensure_equilibrator()
    cc.p_h = Q_(pH)
    cc.ionic_strength = Q_(f"{ionic_strength_M} M")

    results: dict[str, float | None] = {}
    for compound_id, structures in hash_of_structures.items():
        try:
            smiles = _coerce_to_smiles(structures, compound_id)
            compound = _resolve_compound(cc, smiles, compound_id)
            if compound is None:
                print(
                    f"ΔG skip {compound_id}: not in eQuilibrator cache "
                    f"(skeleton not found)"
                )
                results[compound_id] = None
                continue
            dg = cc.standard_dg_prime(Reaction({compound: 1}))
            mean = dg.value.m_as("kJ/mol") if hasattr(dg, "value") else dg.m_as("kJ/mol")
            results[compound_id] = float(mean)
        except Exception as exc:
            print(f"Error calculating ΔG for {compound_id}: {exc}")
            results[compound_id] = None
    return results


def calculate_delta_g_single(
    smiles: str,
    pH: float = 7.0,
    ionic_strength_M: float = 0.25,
) -> tuple[float | None, str | None]:
    """Compute ΔG'° for a single SMILES.

    Returns ``(value_kjmol, error_string_or_None)``.
    """
    cc, Q_, Reaction = _ensure_equilibrator()
    cc.p_h = Q_(pH)
    cc.ionic_strength = Q_(f"{ionic_strength_M} M")
    try:
        if not smiles:
            return None, "empty_smiles"
        compound = _resolve_compound(cc, smiles, "input")
        if compound is None:
            return None, "no_compound_in_cache"
        dg = cc.standard_dg_prime(Reaction({compound: 1}))
        mean = dg.value.m_as("kJ/mol") if hasattr(dg, "value") else dg.m_as("kJ/mol")
        return float(mean), None
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"
