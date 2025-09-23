from typing import Iterable
from enum import Enum
import hashlib
import pathlib
import polars as pl
import numpy as np

def hash_path(rxns: list[tuple[int, str]]) -> str:
    '''
    Returns hash id for pathway given
    reaction ids
    
    Args
    -----
    rxns: list[tuple[int, str]]
        List of (generation, reaction id) tuples
    '''
    rxns = sorted(rxns, key=lambda x: (x[0], x[1])) # Sort by generation, then lexicographically
    reaction_ids = [r[1] for r in rxns]
    concat = "".join(reaction_ids)
    return hashlib.sha1(concat.encode("utf-8")).hexdigest()

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
    
    def __init__(self, study: pathlib.Path, known: pathlib.Path) -> None:
        self.study = study
        self.predicted_reactions = study / "predicted_reactions.parquet"
        self.paths = study / "paths.parquet"
        self.path_stats = study / "path_stats.parquet"
        self.known_reactions = known / "known_reactions.parquet"
        self.enzymes = known / "known_enzymes.parquet"
        self.compounds = study / "compounds.parquet"

        sts = pl.scan_parquet(self.path_stats).select(
            pl.col("starter_ids"),
            pl.col("target_ids"),
            pl.col("starters"),
            pl.col("targets")
        ).explode(
            ["starters", "starter_ids"]
        ).explode(
            ["targets", "target_ids"]
        ).unique()

        starters = sts.select(
            pl.col("starter_ids"),
            pl.col("starters")
        ).unique().collect()

        targets = sts.select(
            pl.col("target_ids"),
            pl.col("targets")
        ).unique().collect()

        self.starters = tuple(starters["starters"])
        self.targets = tuple(targets["targets"])
        self.starter_ids = tuple(starters["starter_ids"])
        self.target_ids = tuple(targets["target_ids"])
        self.starter_name_to_id = dict(zip(starters["starters"], starters["starter_ids"]))
        self.target_name_to_id = dict(zip(targets["targets"], targets["target_ids"]))

    def get_paths(
            self,
            starters: Iterable[str] = None,
            targets: Iterable[str] = None,
            sort_by: str | None = None,
            descending: bool = True,
            lower_bounds: dict[str, float] | None = None,
            upper_bounds: dict[str, float] | None = None,
            filter_by_enzymes: dict[str, Iterable[str]] | None = None,
            top_k: int = None
        ) -> dict[str, pl.DataFrame]:
        '''
        Get valid paths with options to filter & sort

        Args
        -------
        starters: Iterable[str]
            Names of desired starter molecules
        targets: Iterable[str]
            Names of desired target molecules
        sort_by: str
            Field to sort by, can be one of:
                - "mdf"
                - "mean_max_rxn_sim"
                - "mean_mean_rxn_sim"
                - "min_max_rxn_sim"
                - "min_mean_rxn_sim"
                - "feasibility_frac"
        descending: bool = True
            If true sorts in descending order, otherwise ascending
        lower_bounds: dict[str, float]
            Maps field name -> lower bound value.
            Same fields as in sort_by.
        upper_bounds: dict[str, float]
            Maps field name -> upper bound value.
            Same fields as in sort_by.
        filter_by_enzymes: dict[str, Iterable]
            Maps field name -> list of desired values.
            Field may be one of:
                - "id"
                - "sequence"
                - "existence"
                - "reviewed"
                - "ec"
                - "organism"
                - "name"
        top_k: int
            If provided, returns top k paths sorted by sort_by field
        
        Returns
        ---------
        batch: dict[str, pl.DataFrame]
            {
                "paths": pl.DataFrame of paths,
                "predicted_reactions": pl.DataFrame of predicted reactions,
                "known_reactions": pl.DataFrame of known reactions,
                "enzymes": pl.DataFrame of enzymes
            }
        '''
        if top_k and not sort_by:
            raise ValueError("If top_k is provided, sort_by must also be provided")
        
        _path_stats_schema = pl.read_parquet_schema(self.path_stats)
        _enzyme_schema = pl.read_parquet_schema(self.enzymes)

        if sort_by and sort_by not in _path_stats_schema:
            raise ValueError(f"Invalid sort_by field: {sort_by}. Must be one of {_path_stats_schema}")
        
        if lower_bounds and not all([f in _path_stats_schema for f in lower_bounds]):
            raise ValueError(f"Invalid lower_bounds field(s): {lower_bounds}. Must be one of {_path_stats_schema}")
        
        if upper_bounds and not all([f in _path_stats_schema for f in upper_bounds]):
            raise ValueError(f"Invalid upper_bounds field(s): {upper_bounds}. Must be one of {_path_stats_schema}")
        
        if filter_by_enzymes and not all([f in _enzyme_schema for f in filter_by_enzymes]):
            raise ValueError(f"Invalid filter_by_enzymes field(s): {filter_by_enzymes}. Must be one of {_enzyme_schema}")
        
        path_stats_lf = pl.scan_parquet(self.path_stats)
        
        if starters:
            starter_ids = []
            for s in starters:
                if s not in self.starter_name_to_id:
                    raise ValueError(f"Invalid starter name: {s}. Must be one of {self.starters}")
                starter_ids.append(self.starter_name_to_id[s])

            path_stats_lf = path_stats_lf.filter(
                pl.col("starter_ids").list.eval(pl.element().is_in(starter_ids)).list.all()
            )
        
        if targets:
            target_ids = []
            for t in targets:
                if t not in self.target_name_to_id:
                    raise ValueError(f"Invalid target name: {t}. Must be one of {self.targets}")
                target_ids.append(self.target_name_to_id[t])

            path_stats_lf = path_stats_lf.filter(
                pl.col("target_ids").list.eval(pl.element().is_in(target_ids)).list.all()
            )

        if lower_bounds:
            for f, v in lower_bounds.items():
                path_stats_lf = path_stats_lf.filter(pl.col(f) >= v)
        
        if upper_bounds:
            for f, v in upper_bounds.items():
                path_stats_lf = path_stats_lf.filter(pl.col(f) <= v)

        paths_lf = path_stats_lf.join(
            other=pl.scan_parquet(self.paths),
            left_on="id",
            right_on="path_id",
            how="inner"
        )

        if sort_by:
            paths_lf = paths_lf.sort(sort_by, descending=descending)

        paths_df = paths_lf.slice(0, top_k).collect()

        prids = set(paths_df['rxn_id']) # Get all unique reaction ids in paths
        prxns_df = pl.scan_parquet(self.predicted_reactions).filter(pl.col("id").is_in(prids)).collect()
        prxns_df = prxns_df.with_columns(
            pl.col("id").map_elements(lambda x: str(self.study / "svgs" / f"{x}.svg"), return_dtype=pl.String).alias("image"),
        )

        krids = set(prxns_df['analogue_ids'].explode()) # Get all unique known reaction ids in predicted reactions
        krxns_df = pl.scan_parquet(self.known_reactions).filter(pl.col("id").is_in(krids)).collect()
        krxns_df = krxns_df.with_columns(
            pl.col("id").map_elements(lambda x: str(self.study / "svgs" / f"{x}.svg"), return_dtype=pl.String).alias("image"),
        )

        enz_ids = set(krxns_df['enzymes'].explode()) # Get all unique enzyme ids in known reactions
        enz_lf = pl.scan_parquet(self.enzymes).filter(pl.col("id").is_in(enz_ids))

        if filter_by_enzymes:
            for f, v in filter_by_enzymes.items():
                enz_lf = enz_lf.filter(pl.col(f).is_in(v))

        enz_df = enz_lf.collect()

        return {"paths": paths_df, "predicted_reactions": prxns_df, "known_reactions": krxns_df, "enzymes": enz_df}

    # TODO: Reimplement w/ polars
    def get_path_with_id(self, pid: str) -> dict[str, pl.DataFrame]:
        pass

def pick_constraints_for_MDF(S: any, Nc: any, Nr: any, ln_conc: any, dg_prime: any, B: any, lb: float, ub: float) -> list:
    '''
    AUTHOR: YASH CHAINANI
    '''
    
    ### Cofactors to track for specific constraints
    ATP_inchi_key = "ZKHQWZAMYRWXGA-UHFFFAOYSA-N"
    ADP_inchi_key = "XTWYTFMLZFPYCI-UHFFFAOYSA-N"
    AMP_inchi_key = "UDMBCSSLTHHNCD-UHFFFAOYSA-N"
    NADP_plus_inchi_key = "XJLXINKUBYWONI-UHFFFAOYSA-O"
    NADPH_inchi_key = "ACFIXJIJDZMPPO-UHFFFAOYSA-N"
    NAD_plus_inchi_key = "BAWFJGJZGIEFAR-UHFFFAOYSA-O"
    NADH_inchi_key = "BOPGDPNILDQYTO-UHFFFAOYSA-N"
    PI_inchi_key = "NBIIXXVUZAFLBC-UHFFFAOYSA-L"
    PPI_inchi_key = "XPPKVPWEQAFLFU-UHFFFAOYSA-K"
    CoA_inchi_key = "RGJOEKWQDUBAIZ-UHFFFAOYSA-N"
    NH3_inchi_key = "QGZKDVFQNNGYKY-UHFFFAOYSA-O"

    ### Store the inchi keys of all the compounds in the reaction/ pathway
    compound_inchi_keys_list = []
    for compound in list(S.index):
        compound_inchi_keys_list.append(compound.inchi_key)

    ### Note the various conditions

    # Rules 14, 15, 223, 224, 408, 409
    ATP_ADP_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 49, 50, 356, 402,
    ATP_ADP_PI_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
        PI_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 170, 171
    ATP_ADP_PI_CoA_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
        PI_inchi_key in compound_inchi_keys_list,
        CoA_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        AMP_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 66, 67, 345, 346
    ATP_AMP_PPI_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_PPI_CoA_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,  # Rules 38, 39
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
        CoA_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_PPI_NH3_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,  # Rules 385, 386
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
        NH3_inchi_key in compound_inchi_keys_list,
    ]

    NADP_NADPH_cofactors = [
        NADP_plus_inchi_key in compound_inchi_keys_list,
        NADPH_inchi_key in compound_inchi_keys_list,
    ]

    NAD_NADH_cofactors = [
        NAD_plus_inchi_key in compound_inchi_keys_list,
        NADH_inchi_key in compound_inchi_keys_list,
    ]

    ### ATP, ADP: phosphate donor/ acceptor rxns (ATP/ ADP = 10)
    # Rules 14, 15
    if all(ATP_ADP_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == ATP_inchi_key:
                ATP_index = i

            if compound_inchi == ADP_inchi_key:
                ADP_index = i

        ATP_ADP_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[ATP_index] - ln_conc[ADP_index]) <= 2.303,
            2.300 <= (ln_conc[ATP_index] - ln_conc[ADP_index]),
        ]

        return ATP_ADP_constraints

    ### ATP, AMP, PPI, CoA: pyrophosphate donor/ acceptor rxns (ATP/ AMP = 10, [pyrophosphate] = 1 mM)
    # Rules 66, 67, 345, 346
    if all(ATP_AMP_PPI_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == ATP_inchi_key:
                ATP_index = i

            if compound_inchi == AMP_inchi_key:
                AMP_index = i

            if compound_inchi == PPI_inchi_key:
                PPI_index = i

        ATP_AMP_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[ATP_index] - ln_conc[AMP_index]) <= 2.303,
            2.303 <= (ln_conc[ATP_index] - ln_conc[AMP_index]),
            ln_conc[PPI_index] <= -4.6,
            -4.7 <= ln_conc[PPI_index],
        ]

        return ATP_AMP_constraints

    ### NADP+ - NADPH reactions (NADPH/ NADP+ = 10)
    # Rules 2,3..
    if all(NADP_NADPH_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == NADP_plus_inchi_key:
                NADP_plus_index = i

            if compound_inchi == NADPH_inchi_key:
                NADPH_index = i

        NADP_NADPH_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[NADPH_index] - ln_conc[NADP_plus_index]) <= 2.303,
            2.300 <= (ln_conc[NADPH_index] - ln_conc[NADP_plus_index]),
        ]

        return NADP_NADPH_constraints

    ### NAD+ - NADH reactions (NADH/ NAD+ = 0.1, i.e NAD+/ NADH = 10):
    # Rules 2,3,...
    if all(NAD_NADH_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if NAD_plus_inchi_key == compound_inchi:
                NAD_plus_index = i

            if NADH_inchi_key == compound_inchi:
                NADH_index = i

        NAD_NADH_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[NAD_plus_index] - ln_conc[NADH_index]) <= 2.303,
            2.300 <= (ln_conc[NAD_plus_index] - ln_conc[NADH_index]),
        ]
        return NAD_NADH_constraints

    ### All others
    else:
        reg_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
        ]

    return reg_constraints
    
if __name__ == '__main__':
    from hydra import compose, initialize
    
    with initialize(version_base=None, config_path="conf/filepaths"):
        filepaths = compose(config_name="filepaths")