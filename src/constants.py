from enum import StrEnum, unique, Enum

from rdkit import Chem


class MolMixin:
    def as_mol(self):
        return Chem.MolFromSmarts(self.value)


@unique
class BaseSMARTS(MolMixin, StrEnum):
    CARBON = '[#6]'
    CARBON_WITH_FUNC_GROUP_ANY = '[#6;$([#6]=[!#6]),$([#6]#[!#6]),$([#6]-[$(N),$(O),$(S),$(P),$(F),$(Cl),$(Br),$(I)])]'
    RING_BOND = '@'
    NON_RING_BOND = '!@'


@unique
class MoleculeSMARTS(MolMixin, StrEnum):
    C4_FG123 = BaseSMARTS.NON_RING_BOND.join([
        BaseSMARTS.CARBON,
        BaseSMARTS.CARBON_WITH_FUNC_GROUP_ANY,
        BaseSMARTS.CARBON_WITH_FUNC_GROUP_ANY,
        BaseSMARTS.CARBON_WITH_FUNC_GROUP_ANY,
    ])
