from src.pickaxe_processing import pk_rhash_to_smarts
from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict, Iterable
from enum import Enum

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
    
    def filter_enzymes(self, filter_by:Dict[str, Iterable]):
        tmp = []
        for e in self.enzymes:
            include = True
            for k, v in filter_by.items():
                if getattr(e, k) not in v:
                    include = False
                    break
            
            if include:
                tmp.append(e)
        
        self.enzymes = tmp 

    @property
    def has_valid_enzymes(self):
        return True if len(self.enzymes) > 0 else False

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
    
    def filter_analogues(self, filter_by:Dict[str, Iterable]):
        tmp = {}
        for krid, a in self.analogues.items():
            include = True
            for k, v in filter_by.items():
                if getattr(a, k) not in v:
                    include = False
                    break
            
            if include:
                tmp[krid] = a
        
        self.analogues = tmp
    
    @property
    def has_valid_analogues(self):
        return True if len(self.analogues) > 0 else False
    
    def top_rcmcs(self, k:int):
        '''Top k RCMCS scores'''
        srt_krids = sorted(self.rcmcs, key=lambda x : self.rcmcs[x], reverse=True)[:k]
        return [self.rcmcs[krid] for krid in srt_krids]
    
    def top_analogues(self, k:int):
        '''Top k analogues, scored on RCMCS'''
        srt_krids = sorted(self.rcmcs, key=lambda x : self.rcmcs[x], reverse=True)[:k]
        return [self.analogues[krid] for krid in srt_krids]

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
    
    @property
    def valid(self):
        return all([r.has_valid_analogues for r in self.reactions])
    
    @property
    def min_rcmcs(self):
        return min([r.top_rcmcs(k=1)[0] for r in self.reactions])
    
    @property
    def mean_rcmcs(self):
        top_rcmcs = [r.top_rcmcs(k=1)[0] for r in self.reactions]
        return sum(top_rcmcs) / len(top_rcmcs)
    
    def aggregate_rcmcs(self, pr_agg:str, kr_agg:str, k:int):
        aggs = {
            'min': lambda x : min(x),
            'mean': lambda x: sum(x) / len(x),
            'max': lambda x: max(x)
        }
        if pr_agg not in aggs or kr_agg not in aggs:
            raise ValueError(f"Choose valid aggregation methods from: {aggs.keys()}")
        
        top_k_rcmcs = [r.top_rcmcs(k=k) for r in self.reactions]
        kr_agged_rcmcs = [aggs[kr_agg](top) for top in top_k_rcmcs]
        return aggs[pr_agg](kr_agged_rcmcs)        
    
class EnzymeExistence(Enum):
    PROTEIN = 'Evidence at protein level'
    TRANSCRIPT = 'Evidence at transcript level'
    HOMOLOGY = 'Inferred from homology'
    PREDICTED = 'Predicted'
    UNCERTAIN = 'Uncertain'

class PathWrangler:
    enzyme_existence = EnzymeExistence
    def __init__(self) -> None:
        pass

    def retrieve(self):
        pass

    def load(self, path_filepath, pr_filepath, kr_filepath):
        pass


if __name__ == '__main__':
    e1 = Enzyme('up1', 'AARTGW', 'Evidence at protein level', 'reviewed', '1.1.1.1', 'mouse', 'mouse protein')
    e2 = Enzyme('up2', 'QWYPOI', 'Evidence at transcript level', 'reviewed', '1.1.1.2', 'e. coli', 'bacteria protein')
    db1 = DatabaseEntry('rhea', '123')
    db2 = DatabaseEntry('rhea', '456')
    es = [e1, e2]
    dbs = [db1, db2]
    kr = KnownReaction('1', 'CC=O>>CCO', '/home/stef/bottle/artifacts/imgs/img1.svg', es, dbs)
    kr2 = KnownReaction('2', 'CC(=O)O>>CC.O=C=O', '/home/stef/bottle/artifacts/imgs/img2.svg', es, dbs)
    pr = PredictedReaction(
        id='1',
        smarts='O=CCC(=O)O>>O=CCC.O=C=O',
        image='/home/stef/bottle/artifacts/imgs/img3.svg',
        operators=['rule1', 'rule2'],
        analogues={'1':kr, '2': kr2},
        rcmcs={'1': 0.9, '2': 0.2}
    )

    # to / from dict
    kr_dict = asdict(kr)
    kr_from = kr.from_dict(kr_dict)
    pr_dict = pr.to_dict()
    pr_from = PredictedReaction.from_dict(pr_dict)

    # Filter enzymes in known reactions
    kr.filter_enzymes({'existence': [EnzymeExistence.PROTEIN.value]})
    print(kr.has_valid_enzymes)
    kr2.filter_enzymes({'existence': [EnzymeExistence.HOMOLOGY.value]})

    # Filter known reactions in predicted reaction
    pr.filter_analogues({'has_valid_enzymes': [True]})

    # Sort analogues
    pr = PredictedReaction(
        id='1',
        smarts='O=CCC(=O)O>>O=CCC.O=C=O',
        image='/home/stef/bottle/artifacts/imgs/img3.svg',
        operators=['rule1', 'rule2'],
        analogues={'1':kr, '2': kr2},
        rcmcs={'1': 0.9, '2': 0.2}
    )

    print(pr.top_analogues(k=2))
    print(pr.top_analogues(k=1))

    # Top RCMCS
    print(pr.top_rcmcs(k=2))
    print(pr.top_rcmcs(k=1))

    # Path-level operations
    path = Path(
        id='1',
        starter='akg',
        target='hopa',
        reactions=[pr, pr],
        _sid='Csdlfk',
        _tid="Caweor34np"
    )

    print(path.valid)
    print(path.min_rcmcs)
    print(path.mean_rcmcs)
    print(path.aggregate_rcmcs(pr_agg='min', kr_agg='mean', k=5))
    

    print('hold')
