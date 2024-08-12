from src.utils import load_json
from src.pickaxe_processing import pk_rhash_to_smarts
from dataclasses import dataclass, asdict, field
from typing import Optional, List, Dict, Iterable
from enum import Enum
import hashlib

def get_path_id(reaction_ids:Iterable[str]):
    '''Returns hash id for pathway given
    reaction ids'''
    concat = "".join(reaction_ids)
    return "P" + hashlib.sha1(concat.encode("utf-8")).hexdigest()

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
    reaction_center:Iterable[Iterable] = field(default_factory=tuple)
    image:str = ''

    def to_dict(self):
        return asdict(self)
    
    @classmethod
    def from_dict(cls, kr:dict):
        '''
        Construct from dict

        Args
        ----
        kr:dict
            cls.to_dict() output
        '''
        objectifiers = {
            'enzymes': lambda L : [Enzyme.from_dict(e) for e in L],
            'db_entries': lambda L : [DatabaseEntry.from_dict(d) for d in L]
        }
        kwargs = dict((k, objectifiers[k](v)) if k in objectifiers else (k, v) for k, v in kr.items())
        return cls(**kwargs)
    
    def filter_enzymes(self, filter_by:Dict[str, Iterable]):
        '''
        Removes enzymes without desired value in given field

        Args
        ----
        filter_by:Dict[str, Iterable]
            Maps field name -> list of valid values
        '''
        
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
        '''Returns dict representation w/ all fields
        except for analogues, only ids are kept'''
        self_dict = asdict(self)
        self_dict['analogues'] = list(self_dict['analogues'].keys())
        return self_dict

    @classmethod
    def from_dict(cls, pr:dict, stored_known_reactions:Dict[str, KnownReaction]):
        '''
        Construct from stored predicted reaction and stored known reactions

        Args
        -----
        pr:dict
            cls.to_dict() output
        stored_known_reactions:dict
            Contains KnownReaction.to_dict() outputs stored under known 
            reaction ids as keys
        '''
        objectifiers = {
            'analogues': lambda L : {krid : stored_known_reactions[krid] for krid in L} if L else {}
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
        if len(self.analogues) == 0:
            return False
        elif not all([a.has_valid_enzymes for a in self.analogues.values()]):
            return False
        else:
            return True
    
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
        '''Returns dict representation w/ all fields
        except for reactions, only ids are kept'''
        self_dict = asdict(self)
        self_dict['reactions'] = [r.id for r in self.reactions]
        return self_dict

    @classmethod
    def from_dict(cls, path:dict, prs:Dict[str, PredictedReaction]):
        '''
        Construct from stored path, stored predicted reactions, &
        stored known reactions

        Args
        -----
        path:dict
            cls.to_dict() output
        stored_known_reactions:dict
            Contains PredictedReaction.to_dict() outputs stored under predicted 
            reaction ids as keys
        stored_known_reactions:dict
            Contains KnownReaction.to_dict() outputs stored under known 
            reaction ids as keys
        '''
        objectifiers = {
            'reactions': lambda L : [prs[prid] for prid in L] if L else []
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
    '''
    Loads dict representations of stored paths, 
    predicted reactions, & known reactions and provides
    method to get valid paths with options to filter & sort
    based on various criteria. 
    '''
    enzyme_existence = EnzymeExistence
    
    def __init__(self, path_filepath:str, pr_filepath:str, kr_filepath:str) -> None:
        self.known_reactions = load_json(kr_filepath)
        self.predicted_reactions = load_json(pr_filepath)
        self.paths = load_json(path_filepath)
        self.starter_targets = self._extract_starter_targets()
        
    def _extract_starter_targets(self):
        tmp = set()
        for k, v in self.paths.items():
            tmp.add((v['starter'], v['target']))
        return tuple(tmp)

    def get_paths(
            self,
            starters:Iterable[str] = None,
            targets:Iterable[str] = None,
            filter_by_enzymes:Dict[str, Iterable] = None,
            sort_by= None,
        ):
        '''
        Get valid paths with options to filter & sort

        Args
        -------
        starters:Iterable[str]
            Names of desired starter molecules
        targets:Iterable[str]
            Names of desired target molecules
        filter_by_enzymes:Dict[str, Iterable]
            Maps field name -> list of valid values
        sort_by:str | dict
            'min_rcmcs' sorts by path's minimum reaction center mcs
                given by taking max RCMCS for each predicted reaction 
            'mean_rcmcs' sorts by paths' mean reaction center mcs
                given by taking max RCMCS for each predicted reaction
            You may provide a dict that specifies how to aggregate RCMCS on
                both the path level (multiple predicted reactions) and the 
                predicted reaction level (multiple analogues)
                In {'kr_agg': str, 'pr_agg': str, 'k': int}
                pr_agg aggregates on path level 'min' | 'mean' | 'max'
                kr_agg aggregates on predicted reaction level 'min' | 'mean' | 'max'
                k takes k top analogues per predicted reaction

        Returns
        ---------
        Fitlered & sorted List[Path]
        '''
        if starters and targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['starter'] in starters and v['target'] in targets]
        elif starters and not targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['starter'] in starters]
        elif not starters and targets:
            req_path_ids = [v['id'] for v in self.paths.values() if v['target'] in targets]
        else:
            req_path_ids = [v['id'] for v in self.paths.values()]

        if filter_by_enzymes:
            if 'existence' in filter_by_enzymes:
                filter_by_enzymes['existence'] = [self.enzyme_existence[v.upper()].value for v in filter_by_enzymes['existence']]
        
        req_prids = set()
        req_krids = set()
        for pid in req_path_ids:
            for prid in self.paths[pid]['reactions']:
                req_prids.add(prid)
                for krid in self.predicted_reactions[prid]['analogues']:
                    req_krids.add(krid)

        krs = {}
        for krid in req_krids:
            kr = KnownReaction.from_dict(self.known_reactions[krid])
            if filter_by_enzymes:
                kr.filter_enzymes(filter_by=filter_by_enzymes)
            krs[krid] = kr

        prs = {}
        for prid in req_prids:
            pr = PredictedReaction.from_dict(self.predicted_reactions[prid], krs)
            prs[prid] = pr

        paths = []
        for pid in req_path_ids:
            p = Path.from_dict(self.paths[pid], prs)
            if p.valid:
                paths.append(p)

        if sort_by is None:
            return paths
        elif type(sort_by) is str:
            return sorted(paths, key= lambda p : getattr(p, sort_by), reverse=True)
        elif type(sort_by) is dict:
            self._validate_custom_agg(sort_by)
            return sorted(paths, key= lambda p : p.aggregate_rcmcs(**sort_by))
        
    def _validate_custom_agg(self, sort_by:dict):
        req_kwargs = {'kr_agg': str, 'pr_agg': str, 'k': int}
        for k, t in req_kwargs.items():
            if k not in sort_by:
                raise ValueError(f"Missing required argument {k}")
            if type(sort_by[k]) is not t:
                raise ValueError(f"Invalid type {sort_by[k]} for argument {k}. Must be {t}")

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
    pr_from = PredictedReaction.from_dict(pr_dict, {'1':kr, '2':kr2})

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

    # Path wrangling
    path_filepath = '../artifacts/processed_expansions/found_paths.json'
    predicted_reactions_filepath = "../artifacts/processed_expansions/predicted_reactions.json"
    known_reactions_filepath = "../artifacts/processed_expansions/known_reactions.json"
    pw = PathWrangler(
        path_filepath=path_filepath,
        pr_filepath=predicted_reactions_filepath,
        kr_filepath=known_reactions_filepath
    )

    pw.get_paths()
    pw.get_paths(sort_by='min_rcmcs')
    pw.get_paths(sort_by='mean_rcmcs')
    pw.get_paths(filter_by_enzymes={'existence':['protein']})
    
    print('hold')
