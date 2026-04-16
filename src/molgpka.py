from pathlib import Path

import math
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from copy import deepcopy
from functools import lru_cache
from itertools import combinations
from torch.nn import Linear, BatchNorm1d, Parameter
from torch_geometric.nn import GlobalAttention
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.utils import scatter, add_remaining_self_loops
from torch_geometric.data import Data
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


# ---------------------------------------------------------------------------
# Module-level paths (robust regardless of working directory)
# ---------------------------------------------------------------------------

_MODULE_DIR = Path(__file__).parent
SMARTS_FILE = _MODULE_DIR / "data" / "smarts_pattern.tsv"
BASE_MODEL_PATH = _MODULE_DIR / "models" / "weight_base.pth"
ACID_MODEL_PATH = _MODULE_DIR / "models" / "weight_acid.pth"
# Optional ModelSeed-trained weights from OPAM2; not required for the stock
# MolGpKa workflow and are only used by callers that opt into them.
MODELSEED_BASE_MODEL_PATH = _MODULE_DIR / "models" / "weight_base_modelseed.pth"
MODELSEED_ACID_MODEL_PATH = _MODULE_DIR / "models" / "weight_acid_modelseed.pth"

# Soft cap used by ``cap_unstable_sites`` when callers want to bound
# combinatorial protonation enumeration: 2**MAX_UNSTABLE_SITES states max.
MAX_UNSTABLE_SITES = 8

n_features = 29
hidden = 1024


# ---------------------------------------------------------------------------
# GCN architecture
# ---------------------------------------------------------------------------

def glorot(tensor: torch.Tensor | None) -> None:
    if tensor is not None:
        stdv = math.sqrt(6.0 / (tensor.size(-2) + tensor.size(-1)))
        tensor.data.uniform_(-stdv, stdv)


def zeros(tensor: torch.Tensor | None) -> None:
    if tensor is not None:
        tensor.data.fill_(0)


class GCNConv(MessagePassing):
    r"""Graph convolutional operator from `"Semi-supervised Classification
    with Graph Convolutional Networks" <https://arxiv.org/abs/1609.02907>`_.

    .. math::
        \mathbf{X}^{\prime} = \mathbf{\hat{D}}^{-1/2} \mathbf{\hat{A}}
        \mathbf{\hat{D}}^{-1/2} \mathbf{X} \mathbf{\Theta},

    where :math:`\mathbf{\hat{A}} = \mathbf{A} + \mathbf{I}` denotes the
    adjacency matrix with inserted self-loops and
    :math:`\hat{D}_{ii} = \sum_{j=0} \hat{A}_{ij}` its diagonal degree matrix.

    Args:
        in_channels: Size of each input sample.
        out_channels: Size of each output sample.
        improved: If True, computes :math:`\mathbf{\hat{A}}` as
            :math:`\mathbf{A} + 2\mathbf{I}`. Default: False.
        cached: If True, caches the normalised adjacency on first forward pass.
            Only suitable for transductive learning. Default: False.
        bias: If False, no additive bias. Default: True.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        improved: bool = False,
        cached: bool = False,
        bias: bool = True,
        **kwargs,
    ) -> None:
        super(GCNConv, self).__init__(aggr="add", **kwargs)

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.improved = improved
        self.cached = cached

        self.weight = Parameter(torch.Tensor(in_channels, out_channels))

        if bias:
            self.bias = Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter("bias", None)

        self.reset_parameters()

    def reset_parameters(self) -> None:
        glorot(self.weight)
        zeros(self.bias)
        self.cached_result = None
        self.cached_num_edges = None

    @staticmethod
    def norm(edge_index, num_nodes, edge_weight=None, improved=False, dtype=None):
        if edge_weight is None:
            edge_weight = torch.ones(
                (edge_index.size(1),), dtype=dtype, device=edge_index.device
            )

        fill_value = 1 if not improved else 2
        edge_index, edge_weight = add_remaining_self_loops(
            edge_index, edge_weight, fill_value, num_nodes
        )

        row, col = edge_index
        deg = scatter(edge_weight, row, dim=0, dim_size=num_nodes, reduce="sum")
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float("inf")] = 0

        return edge_index, deg_inv_sqrt[row] * edge_weight * deg_inv_sqrt[col]

    def forward(self, x, edge_index, edge_weight=None):
        x = torch.matmul(x, self.weight)

        if self.cached and self.cached_result is not None:
            if edge_index.size(1) != self.cached_num_edges:
                raise RuntimeError(
                    "Cached {} number of edges, but found {}. Please "
                    "disable the caching behavior of this layer by removing "
                    "the `cached=True` argument in its constructor.".format(
                        self.cached_num_edges, edge_index.size(1)
                    )
                )

        if not self.cached or self.cached_result is None:
            self.cached_num_edges = edge_index.size(1)
            edge_index, norm = self.norm(
                edge_index, x.size(0), edge_weight, self.improved, x.dtype
            )
            self.cached_result = edge_index, norm

        edge_index, norm = self.cached_result
        return self.propagate(edge_index, x=x, norm=norm)

    def message(self, x_j, norm):
        return norm.view(-1, 1) * x_j

    def update(self, aggr_out):
        if self.bias is not None:
            aggr_out = aggr_out + self.bias
        return aggr_out

    def __repr__(self) -> str:
        return "{}({}, {})".format(
            self.__class__.__name__, self.in_channels, self.out_channels
        )


class GCNNet(torch.nn.Module):
    """5-layer GCN for per-atom pKa prediction with global attention pooling."""

    def __init__(self) -> None:
        super(GCNNet, self).__init__()
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

        self.att = GlobalAttention(Linear(hidden, 1))
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
        x = F.relu(self.conv1(x, edge_index))
        x = self.bn1(x)
        x = F.relu(self.conv2(x, edge_index))
        x = self.bn2(x)
        x = F.relu(self.conv3(x, edge_index))
        x = self.bn3(x)
        x = F.relu(self.conv4(x, edge_index))
        x = self.bn4(x)
        x = F.relu(self.conv5(x, edge_index))
        x = self.bn5(x)
        x = self.att(x, batch)
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        return self.fc4(x)


# ---------------------------------------------------------------------------
# Molecule featurisation
# ---------------------------------------------------------------------------

def get_bond_pair(mol: Chem.Mol) -> list[list[int]]:
    """Return bidirectional edge index for all bonds in *mol*."""
    bonds = mol.GetBonds()
    res: list[list[int]] = [[], []]
    for bond in bonds:
        res[0] += [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
        res[1] += [bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()]
    return res


def one_hot(x, allowable_set: list) -> list[bool]:
    """One-hot encode *x* against *allowable_set*, mapping unknowns to the last entry."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def get_atom_features(mol: Chem.Mol, aid: int) -> list[list]:
    """Compute a 29-dimensional feature vector for every atom in *mol*.

    The target atom is identified by *aid* (atom index) and receives a
    distance-to-self of 0 and a flag marking it as the ionisation site.
    """
    AllChem.ComputeGasteigerCharges(mol)
    Chem.AssignStereochemistry(mol)

    acceptor_smarts_one = "[!$([#1,#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]"
    acceptor_smarts_two = "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"
    donor_smarts_one = "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]"
    donor_smarts_two = "[!$([#6,H0,-,-2,-3]),$([!H0;#7,#8,#9])]"

    hydrogen_donor_match = list(
        set(
            mol.GetSubstructMatches(Chem.MolFromSmarts(donor_smarts_one))
            + mol.GetSubstructMatches(Chem.MolFromSmarts(donor_smarts_two))
        )
    )
    hydrogen_acceptor_match = list(
        set(
            mol.GetSubstructMatches(Chem.MolFromSmarts(acceptor_smarts_one))
            + mol.GetSubstructMatches(Chem.MolFromSmarts(acceptor_smarts_two))
        )
    )

    ring = mol.GetRingInfo()
    m = []
    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        o: list = []
        o += one_hot(
            atom.GetSymbol(), ["C", "H", "O", "N", "S", "Cl", "F", "Br", "P", "I"]
        )
        o += [atom.GetDegree()]
        o += one_hot(
            atom.GetHybridization(),
            [
                Chem.rdchem.HybridizationType.SP,
                Chem.rdchem.HybridizationType.SP2,
                Chem.rdchem.HybridizationType.SP3,
                Chem.rdchem.HybridizationType.SP3D,
                Chem.rdchem.HybridizationType.SP3D2,
            ],
        )
        o += [atom.GetImplicitValence()]
        o += [atom.GetIsAromatic()]
        o += [
            ring.IsAtomInRingOfSize(atom_idx, 3),
            ring.IsAtomInRingOfSize(atom_idx, 4),
            ring.IsAtomInRingOfSize(atom_idx, 5),
            ring.IsAtomInRingOfSize(atom_idx, 6),
            ring.IsAtomInRingOfSize(atom_idx, 7),
            ring.IsAtomInRingOfSize(atom_idx, 8),
        ]
        o += [atom_idx in hydrogen_donor_match]
        o += [atom_idx in hydrogen_acceptor_match]
        o += [atom.GetFormalCharge()]
        o += [0 if atom_idx == aid else len(Chem.rdmolops.GetShortestPath(mol, atom_idx, aid))]
        o += [atom_idx == aid]
        m.append(o)
    return m


def mol2vec(
    mol: Chem.Mol,
    atom_idx: int,
    evaluation: bool = True,
    pka: float | None = None,
) -> Data:
    """Convert *mol* and ionisation site *atom_idx* into a PyG ``Data`` object."""
    node_f = get_atom_features(mol, atom_idx)
    edge_index = get_bond_pair(mol)
    if evaluation:
        batch = np.zeros(len(node_f))
        return Data(
            x=torch.tensor(node_f, dtype=torch.float32),
            edge_index=torch.tensor(edge_index, dtype=torch.long),
            batch=torch.tensor(batch, dtype=torch.long),
        )
    return Data(
        x=torch.tensor(node_f, dtype=torch.float32),
        edge_index=torch.tensor(edge_index, dtype=torch.long),
        pka=torch.tensor([[pka]], dtype=torch.float),
    )


# ---------------------------------------------------------------------------
# Model loading and inference
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4)
def _load_model_cached(model_file: str, device: str) -> GCNNet:
    """Cached inner loader; keyed on stringified path + device."""
    model = GCNNet().to(device)
    # ``weights_only=True`` is safe for standard state-dicts and guards
    # against malicious pickles; matches the OPAM2 refinement.
    try:
        state = torch.load(model_file, map_location=device, weights_only=True)
    except TypeError:
        # Older PyTorch (< 2.0) does not support weights_only.
        state = torch.load(model_file, map_location=device)
    model.load_state_dict(state)
    model.eval()
    return model


def load_model(model_file: Path | str, device: str = "cpu") -> GCNNet:
    """Load a ``GCNNet`` from a saved state-dict at *model_file*.

    Back-compatible with the original MolGpKa flow: returns a freshly-evaluated
    ``GCNNet``. Underlying deserialisation is cached by ``(path, device)`` so
    repeated calls do not re-read the ``.pth`` file from disk.
    """
    return _load_model_cached(str(model_file), device)


def model_pred(mol: Chem.Mol, aid: int, model: GCNNet, device: str = "cpu") -> float:
    """Run GCN inference and return the predicted pKa for atom *aid* in *mol*."""
    data = mol2vec(mol, aid)
    with torch.no_grad():
        data = data.to(device)
        pka_tensor = model(data)
    return float(pka_tensor.cpu().numpy()[0][0])


# ---------------------------------------------------------------------------
# SMARTS-based ionisable-site detection
# ---------------------------------------------------------------------------

def split_acid_base_pattern(
    smarts_file: Path | str = SMARTS_FILE,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and split the SMARTS pattern file into acid and base DataFrames."""
    df = pd.read_csv(smarts_file, sep="\t")
    return df[df.Acid_or_base == "A"], df[df.Acid_or_base == "B"]


def unique_acid_match(matches: list) -> list:
    """Deduplicate substructure match tuples, keeping single- and double-atom matches."""
    single_matches = [[j] for j in set(m[0] for m in matches if len(m) == 1)]
    double_matches = [m for m in matches if len(m) == 2]
    return double_matches + single_matches


def match_acid(df_smarts_acid: pd.DataFrame, mol: Chem.Mol) -> list[int]:
    """Return atom indices of all acidic sites in *mol*."""
    matches = []
    for _, _name, smarts, index, _ab in df_smarts_acid.itertuples():
        pattern = Chem.MolFromSmarts(smarts)
        match = mol.GetSubstructMatches(pattern)
        if not match:
            continue
        if len(index) > 2:
            indices = [int(i) for i in index.split(",")]
            for m in match:
                matches.append([m[indices[0]], m[indices[1]]])
        else:
            idx = int(index)
            for m in match:
                matches.append([m[idx]])
    flat = []
    for group in unique_acid_match(matches):
        flat.extend(group)
    return flat


def match_base(df_smarts_base: pd.DataFrame, mol: Chem.Mol) -> list[int]:
    """Return atom indices of all basic sites in *mol*."""
    matches = []
    for _, _name, smarts, indexs, _ab in df_smarts_base.itertuples():
        pattern = Chem.MolFromSmarts(smarts)
        match = mol.GetSubstructMatches(pattern)
        if not match:
            continue
        for index in indexs.split(","):
            idx = int(index)
            for m in match:
                matches.append([m[idx]])
    flat = []
    for group in unique_acid_match(matches):
        flat.extend(group)
    return flat


def get_ionization_aid(
    mol: Chem.Mol,
    acid_or_base: str | None = None,
    smarts_file: Path | str = SMARTS_FILE,
) -> list[int] | tuple[list[int], list[int]]:
    """Identify ionisable atom indices in *mol*.

    Parameters
    ----------
    mol:
        RDKit molecule (should have explicit Hs added).
    acid_or_base:
        ``"acid"`` → return only acidic sites; ``"base"`` → only basic;
        ``None`` → return ``(acid_matches, base_matches)``.
    smarts_file:
        Path to the TSV file of SMARTS patterns.
    """
    if mol is None:
        raise RuntimeError("read mol error")
    df_acid, df_base = split_acid_base_pattern(smarts_file)
    acid_matches = match_acid(df_acid, mol)
    base_matches = match_base(df_base, mol)
    if acid_or_base is None:
        return acid_matches, base_matches
    if acid_or_base == "acid":
        return acid_matches
    return base_matches


# ---------------------------------------------------------------------------
# Protonation-state manipulation
# ---------------------------------------------------------------------------

def modify_acid(at: Chem.Atom) -> None:
    """Set atom *at* to its deprotonated (anionic) form."""
    hnum = at.GetNumExplicitHs()
    at.SetFormalCharge(-1)
    at.SetNumExplicitHs(max(hnum - 1, 0))


def modify_base(at: Chem.Atom) -> None:
    """Set atom *at* to its protonated (cationic) form."""
    hnum = at.GetNumExplicitHs()
    at.SetFormalCharge(1)
    at.SetNumExplicitHs(hnum + 1)


def modify_mol(
    mol: Chem.Mol,
    acid_dict: dict[int, float],
    base_dict: dict[int, float],
) -> Chem.Mol:
    """Tag each atom in *mol* with its ionisation type and predicted pKa.

    Returns the molecule with explicit Hs removed.
    """
    for at in mol.GetAtoms():
        idx = at.GetIdx()
        if idx in acid_dict:
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
    """Classify ionisable atoms as stable or unstable relative to *ph* ± *tph*.

    Returns
    -------
    stable_data:
        List of ``[atom_idx, pKa, "A"|"B"]`` for atoms unambiguously in one state.
    unstable_data:
        List of ``[atom_idx, pKa, "A"|"B"]`` for atoms near the pH boundary.
    pkas:
        Dict mapping atom index to pKa value for all ionisable atoms.
    """
    stable_data, unstable_data, pkas = [], [], {}
    for at in mol.GetAtoms():
        props = at.GetPropsAsDict()
        acid_or_basic = props.get("ionization", False)
        if not acid_or_basic or acid_or_basic == "O":
            continue
        pka = float(props.get("pKa", float("nan")))
        idx = at.GetIdx()
        pkas[idx] = pka
        if acid_or_basic == "A":
            if pka < ph - tph:
                stable_data.append([idx, pka, "A"])
            elif ph - tph <= pka <= ph + tph:
                unstable_data.append([idx, pka, "A"])
        elif acid_or_basic == "B":
            if pka > ph + tph:
                stable_data.append([idx, pka, "B"])
            elif ph - tph <= pka <= ph + tph:
                unstable_data.append([idx, pka, "B"])
    return stable_data, unstable_data, pkas


def modify_stable_pka(new_mol: Chem.Mol, stable_data: list) -> None:
    """Apply stable protonation states to *new_mol* in-place."""
    for idx, _pka, acid_or_basic in stable_data:
        at = new_mol.GetAtomWithIdx(idx)
        if acid_or_basic == "A":
            modify_acid(at)
        elif acid_or_basic == "B":
            modify_base(at)


def modify_unstable_pka(
    mol: Chem.Mol,
    unstable_data: list,
    i: int,
) -> list[str]:
    """Return SMILES for all combinations of *i* unstable atoms being modified.

    Parameters
    ----------
    mol:
        Molecule with stable modifications already applied.
    unstable_data:
        List of ``[atom_idx, pKa, "A"|"B"]`` for boundary-zone atoms.
    i:
        Number of unstable atoms to simultaneously modify.
    """
    new_unsmis = []
    for pka_datas in combinations(unstable_data, i):
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
        new_unsmis.append(smi)
    return new_unsmis


# ---------------------------------------------------------------------------
# OPAM2 refinements (additive — original callers unaffected)
# ---------------------------------------------------------------------------

def prepare_mol_for_prediction(mol: Chem.Mol, uncharged: bool = True) -> Chem.Mol:
    """Normalise *mol* for prediction while preserving heavy-atom indices.

    The older ``predict_for_protonate`` flow ran
    ``MolFromSmiles(MolToSmiles(mol))`` after uncharging, which canonicalised
    atom order and broke index-based comparisons against the input.
    ``Uncharger`` alone preserves atom order, and ``AddHs`` appends hydrogens
    without renumbering heavy atoms, so heavy-atom indices in the returned
    mol match those of the input.
    """
    if uncharged:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    return AllChem.AddHs(mol)


def cap_unstable_sites(
    unstable_data: list,
    ph: float,
    max_sites: int = MAX_UNSTABLE_SITES,
) -> tuple[list, list]:
    """Cap combinatorial enumeration by keeping only *max_sites* borderline
    atoms closest to *ph*.

    Overflow sites are promoted to stable states based on which side of *ph*
    their pKa falls on; this keeps ``2**len(unstable)`` bounded for highly
    ionisable molecules.

    Returns
    -------
    kept, promoted_stable:
        ``kept`` is the truncated unstable list; ``promoted_stable`` is the
        list of overflow sites re-categorised as stable.
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
    acid_model: Path | str | None = None,
    base_model: Path | str | None = None,
    output_format: str = "inchi",
    max_unstable_sites: int = MAX_UNSTABLE_SITES,
    uncharged: bool = True,
    device: str = "cpu",
) -> list:
    """Enumerate protonation microstates at ``ph ± tph`` (OPAM2 entry point).

    Accepts a SMILES *or* an InChI string. Neutralises the input before
    pKa prediction so pre-charged zwitterions are handled correctly, then
    re-assigns charges based on predicted pKa vs target pH.

    Reuses the existing :func:`modify_mol`, :func:`get_pKa_data`,
    :func:`modify_stable_pka`, :func:`modify_unstable_pka`,
    :func:`get_ionization_aid`, :func:`model_pred`, and :func:`load_model`
    primitives from this module — only the orchestration is new.

    Caps combinatorial expansion at *max_unstable_sites* borderline sites
    (default 8 → 256 states) via :func:`cap_unstable_sites`.
    """
    omol = Chem.MolFromSmiles(smi_or_inchi)
    if omol is None:
        omol = Chem.MolFromInchi(smi_or_inchi)
    if omol is None:
        raise ValueError(
            f"Could not parse input as SMILES or InChI: {smi_or_inchi!r}"
        )

    omol = prepare_mol_for_prediction(omol, uncharged=uncharged)

    acid_path = acid_model if acid_model is not None else ACID_MODEL_PATH
    base_path = base_model if base_model is not None else BASE_MODEL_PATH
    acid_net = load_model(acid_path, device=device)
    base_net = load_model(base_path, device=device)
    acid_dict = {
        aid: model_pred(omol, aid, acid_net, device=device)
        for aid in get_ionization_aid(omol, acid_or_base="acid")
    }
    base_dict = {
        aid: model_pred(omol, aid, base_net, device=device)
        for aid in get_ionization_aid(omol, acid_or_base="base")
    }

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
