import re
from enum import StrEnum, unique

from rdkit import Chem


class MolMixin:
    def as_mol(self):
        return Chem.MolFromSmarts(self.value)


@unique
class BaseSmarts(MolMixin, StrEnum):
    CARBON_ANY = '[#6]'
    CARBON_NO_FUNC_GROUP = '[#6;!$([#6]~[!#1&!#6]);!$([#6]@[#6])]'
    REACTIVE_NONMETALS_NO_C = '[#7,#8,#9,#15,#16,#17,#34,#35,#53]'
    CARBON_WITH_FUNC_GROUP_NONMETAL = '[#6;$([#6]~[#7,#8,#9,#15,#16,#17,#34,#35,#53])]'
    # CARBON_WITH_FUNC_GROUP_NONMETAL = '[#6;!$([#6]~[#6]),$([#6]~[#7,#8,#9,#15,#16,#17,#34,#35,#53]~[!#6])]'
    RING_BOND = '@'
    NON_RING_BOND = '!@'


def smarts_activation(
        func_groups: list[str],
        *,
        c_bonds: str = BaseSmarts.NON_RING_BOND,
        max_carbon: int = -1
) -> str:
    """Construct a SMARTS pattern for different activation scenarios.

    Args:
        func_groups (list[str]): ordered list of non-func-group and func-group, all ptns must be enclosed in []
        c_bonds (str): types of bonds to consider (defaults to non-ring)
        max_carbon (int): default implying the same upper bound as the member count of func_groups.

    Returns:
        str: smarts pattern
    """
    n_max = max(len(func_groups), max_carbon)
    ptn_c_upper_bound = f'[!$({c_bonds.join([BaseSmarts.CARBON_ANY] * (n_max + 1))});'
    ptn_prelim = c_bonds.join(func_groups)
    return re.sub('^\\[', ptn_c_upper_bound, ptn_prelim)


@unique
class ActivationSmarts(MolMixin, StrEnum):
    C4_FG123 = smarts_activation([
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL,
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL,
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL,
        BaseSmarts.CARBON_NO_FUNC_GROUP,
    ])
    C5_FG125 = smarts_activation([
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL,
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL,
        BaseSmarts.CARBON_NO_FUNC_GROUP,
        BaseSmarts.CARBON_NO_FUNC_GROUP,
        BaseSmarts.CARBON_WITH_FUNC_GROUP_NONMETAL
    ])
