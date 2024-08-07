from src.pickaxe_processing import pk_rhash_to_smarts
from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict, Iterable

@dataclass
class DatabaseEntry:
    name:str
    id:str

    @classmethod
    def from_dict(cls, dbe:dict):
        return cls(**dbe)

@dataclass
class Enzyme:
    uniprot_id:str
    sequence:Optional[str] = None
    existence:Optional[str] = None
    reviewed:Optional[str] = None
    ec:Optional[str] = None
    organism:Optional[str] = None
    name:Optional[str] = None

    def to_dict(self):
        return asdict(self)
    
    @classmethod
    def from_dict(cls, enz:dict):
        return cls(**enz)

@dataclass
class KnownReaction:
    id:str
    smarts:str
    operators:List[str]
    enzymes:List[Enzyme]
    db_entries:List[DatabaseEntry]
    image:str = ''

    def to_dict(self):
        return asdict(self)
    
    @classmethod
    def from_dict(cls, kr:dict):
        objectifiers = {
            'enzymes': lambda L : [Enzyme.from_dict(e) for e in L],
            'db_entries': lambda L : [DatabaseEntry.from_dict(d) for d in L]
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in kr.items())
        return cls(**kwargs)

@dataclass
class PredictedReaction:
    id:str
    smarts:str
    operators:List[str]
    reaction_center:Iterable[Iterable] = field(default_factory=tuple)
    analogues:Dict[str, KnownReaction] = field(default_factory=dict)
    rcmcs:Dict[str, float] = field(default_factory=dict)
    image:str = ''

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, pr:dict):
        objectifiers = {
            'analogues': lambda D : {k : KnownReaction.from_dict(v) for k, v in D.items()} if D else D
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in pr.items())
        return cls(**kwargs)
    
    @classmethod
    def from_pickaxe(cls, pk, id):
        smarts = pk_rhash_to_smarts(id, pk)
        operators = list(pk.reactions[id]["Operators"])
        return cls(id, smarts, operators)

@dataclass
class Path:
    id:str
    starter:str
    target:str
    reactions:List[PredictedReaction]
    _sid:str # Starter hash
    _tid:str # Target hash

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, path:dict):
        objectifiers = {
            'reactions': lambda L : [PredictedReaction.from_dict(pr) for pr in L]
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in path.items())
        return cls(**kwargs)

if __name__ == '__main__':
    e1 = Enzyme('up1', 'AARTGW', 'Protein', 'reviewed', '1.1.1.1', 'mouse', 'mouse protein')
    e2 = Enzyme('up2', 'QWYPOI', 'Protein', 'reviewed', '1.1.1.2', 'e. coli', 'bacteria protein')
    db1 = DatabaseEntry('rhea', '123')
    db2 = DatabaseEntry('rhea', '456')
    es = [e1, e2]
    dbs = [db1, db2]
    kr = KnownReaction('1', 'CC=O>>CCO', '/home/stef/bottle/artifacts/imgs/img1.svg', es, dbs)
    kr_dict = asdict(kr)
    kr_from = kr.from_dict(kr_dict)
    kr2 = KnownReaction('2', 'CC(=O)O>>CC.O=C=O', '/home/stef/bottle/artifacts/imgs/img2.svg', es, dbs)
    pr = PredictedReaction('1', 'O=CCC(=O)O>>O=CCC.O=C=O', '/home/stef/bottle/artifacts/imgs/img3.svg', ['rule1', 'rule2'], {'1':kr, '2': kr2})
    pr_dict = pr.to_dict()
    pr_from = PredictedReaction.from_dict(pr_dict)


    print('hold')
