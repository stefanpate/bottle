import hydra
from omegaconf import DictConfig
from minedatabase.pickaxe import Pickaxe
from pathlib import Path
import polars as pl
from src.schemas import expansion_reactions_schema, compounds_schema, gens_schema
from tqdm import tqdm
from itertools import chain
from src.network import de_am
from ergochemics.mapping import rc_to_nest
from ergochemics.standardize import hash_compound

def reverse_rules(rules: list[str], rule2rxn: dict[int, str], rxn2rule: dict[str, int], reaction_reverses: dict[str, str]) -> list[str]:
    rules = [int(rule) for rule in rules]
    rev_krs = chain(*[rule2rxn.get(rule, []) for rule in rules])
    fwd_krs = [reaction_reverses[rid] for rid in rev_krs if rid in reaction_reverses]
    fwd_rules = list(set([rxn2rule[rid] for rid in fwd_krs if rid in rxn2rule]))
    if len(fwd_rules) == 0:
        # If no forward rules found, mark that these are reverse rules
        # which will not turn up analogues later on but makes clear why not
        fwd_rules = [f"{rule}_rev" for rule in rules]

    return fwd_rules

def parse_reactions(reactions: dict, rule2rxn: dict[int, str], rxn2rule: dict[str, int], rule2size: dict[int, int], reaction_reverses: dict[str, str], mode: str, rule_set: str) -> list[tuple[str, list[str]]]:
    rxn_data = []
    for rxn in tqdm(reactions.values(), desc="Processing reactions", total=len(reactions)):
        if 'am_rxn' not in rxn:
            continue
        
        rules = [int(elt.split('_')[0]) for elt in rxn['Operators']]
        
        if mode == "retro":
            am_smarts = ">>".join(rxn['am_rxn'].split('>>')[::-1])
            rules = reverse_rules(rules, rule2rxn, rxn2rule, reaction_reverses)
        elif mode == "forward":
            am_smarts = rxn['am_rxn']
        else:
            raise ValueError(f"mode must be 'forward' or 'retro', got {mode}")

        lhs, rhs = de_am(am_smarts)
        smarts = f"{".".join(lhs)}>>{".".join(rhs)}"
        size = max([rule2size.get(rule, 0) for rule in rules])
        rules = [f"{rule_set}:{rule}" for rule in rules]
        rxn_data.append((smarts, am_smarts, rules, mode, size))

    return rxn_data

@hydra.main(version_base=None, config_path="../conf", config_name="parse_expansion")
def main(cfg: DictConfig) -> None:
    
    if cfg.fwd_exp is not None and cfg.rev_exp is not None:
        steps = int(cfg.fwd_exp[0]) + int(cfg.rev_exp[0])
        fwd_rules = cfg.fwd_exp.split('_rules_')[1]
        rev_rules = cfg.rev_exp.split('_rules_')[1]
        starter_name = cfg.fwd_exp.split('_steps_')[1].split('_to_')[0]
        target_name = cfg.rev_exp.split('_steps_')[1].split('_to_')[0]
        exp_name = Path(f"{steps}_steps_{starter_name}_to_{target_name}_combo_{fwd_rules}_and_{rev_rules}_rules")
    elif cfg.fwd_exp is not None:
        exp_name = Path(cfg.fwd_exp).stem
    elif cfg.rev_exp is not None:
        exp_name = Path(cfg.rev_exp).stem
  
    if not Path(exp_name).exists():
        Path(exp_name).mkdir(parents=True, exist_ok=True)
    
    kcs = pl.read_parquet(
        Path(cfg.filepaths.known) / "known_compounds.parquet",
    )
    kc_smi2name = dict(zip(kcs["smiles"], kcs["name"].to_list()))

    rxn_data = []
    pred_kcs = {}
    helpers = {}
    sources = {}
    targets = {}
    checkpoints = {}
    fwd_gens = 0
    rev_gens = 0
    if cfg.fwd_exp is not None:
        fwd_rule_set = cfg.fwd_exp.split('_rules_')[1]
        fwd_exp = Pickaxe()
        fwd_exp.load_pickled_pickaxe(
            Path(cfg.filepaths.raw_data) / Path(cfg.fwd_exp)
        )
        fwd_gens = fwd_exp.generation

        mapped_rxns = pl.read_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / f"mapped_known_reactions_x_{fwd_rule_set}_rules.parquet",
        )

        rule2size = mapped_rxns.group_by("rule_id").agg(
            pl.col("template_aidxs").first().map_elements(lambda x: sum(len(elt) for elt in rc_to_nest(x)[0]), return_dtype=pl.Int32).alias("size")
        )

        rule2size = dict(zip(rule2size["rule_id"], rule2size["size"].to_list()))

        for cid, cpd in fwd_exp.compounds.items():
            if cpd['SMILES'] in kc_smi2name and cpd['Type'] != 'Starting Compound':
                _id = hash_compound(cpd['SMILES'])
                pred_kcs[_id] = (cpd['SMILES'], kc_smi2name[cpd['SMILES']])

            if cid[0] == 'X':
                _id = hash_compound(cpd['SMILES'])
                helpers[_id] = (cpd['SMILES'], cpd.get("ID", None))

            if cpd['Type'] == 'Starting Compound':
                _id = hash_compound(cpd['SMILES'])
                sources[_id] = (cpd['SMILES'], cpd.get("ID", None))
        
        for cpd in fwd_exp.targets.values():
            _id = hash_compound(cpd["SMILES"])
            targets[_id] = (cpd["SMILES"], cpd.get("ID", None))

        rxn_data += parse_reactions(
            fwd_exp.reactions,
            rule2rxn={},
            rxn2rule={},
            rule2size=rule2size,
            reaction_reverses={},
            mode="forward",
            rule_set=fwd_rule_set,
        )
  
    if cfg.rev_exp is not None:
        rev_rule_set = cfg.rev_exp.split('_rules_')[1]
        rev_exp = Pickaxe()
        rev_exp.load_pickled_pickaxe(
            Path(cfg.filepaths.raw_data) / Path(cfg.rev_exp)
        )
        rev_gens = rev_exp.generation

        krs = pl.read_parquet(
            Path(cfg.filepaths.known) / "known_reactions.parquet",
        )
        reaction_reverses = dict(zip(krs["id"], krs["reverse"]))

        mapped_rxns = pl.read_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / f"mapped_known_reactions_x_{rev_rule_set}_rules.parquet",
        )

        rule2rxn = mapped_rxns.group_by("rule_id").agg(
            pl.col("rxn_id")
        )

        rule2size = mapped_rxns.group_by("rule_id").agg(
            pl.col("template_aidxs").first().map_elements(lambda x: sum(len(elt) for elt in rc_to_nest(x)[0]), return_dtype=pl.Int32).alias("size")
        )

        rule2size = dict(zip(rule2size["rule_id"], rule2size["size"].to_list()))
        rule2rxn = dict(zip(rule2rxn["rule_id"], rule2rxn["rxn_id"].to_list()))
        rxn2rule = dict(zip(mapped_rxns["rxn_id"], mapped_rxns["rule_id"].to_list()))

        for cid, cpd in rev_exp.compounds.items():
            if cpd['SMILES'] in kc_smi2name and cpd['Type'] != 'Starting Compound':
                _id = hash_compound(cpd['SMILES'])
                pred_kcs[_id] = (cpd['SMILES'], kc_smi2name[cpd['SMILES']])

            if cid[0] == 'X':
                _id = hash_compound(cpd['SMILES'])
                helpers[_id] = (cpd['SMILES'], cpd.get("ID", None))

            if cpd['Type'] == 'Starting Compound': # in reverse exp, these are the targets
                _id = hash_compound(cpd['SMILES'])
                targets[_id] = (cpd['SMILES'], cpd.get("ID", None))

        for cpd in fwd_exp.targets.values(): # in reverse exp, these are the sources
            _id = hash_compound(cpd["SMILES"])
            sources[_id] = (cpd["SMILES"], cpd.get("ID", None))

        rxn_data += parse_reactions(
            rev_exp.reactions,
            rule2rxn=rule2rxn,
            rxn2rule=rxn2rule,
            rule2size=rule2size,
            reaction_reverses=reaction_reverses,
            mode="retro",
            rule_set=rev_rule_set,
        )

    if cfg.fwd_exp is not None and cfg.rev_exp is not None:
        for k in fwd_exp.compounds.keys() & rev_exp.compounds.keys():
            cpd = fwd_exp.compounds[k]
            _id = hash_compound(cpd['SMILES'])
            checkpoints[_id] = (cpd['SMILES'], cpd.get("ID", None))

    # User-defined sources and targets override all other designations
    # helper designation overrides known and checkpoint designations
    checkpoints = {k: v for k, v in checkpoints.items() if k not in sources and k not in targets and k not in helpers}
    pred_kcs = {k: v for k, v in pred_kcs.items() if k not in sources and k not in targets and k not in helpers}
    helpers = {k: v for k, v in helpers.items() if k not in sources and k not in targets}

    cpds = []
    for k, (smiles, name) in sources.items():
        cpds.append((k, smiles, "source", name))
    for k, (smiles, name) in targets.items():
        cpds.append((k, smiles, "target", name))
    for k, (smiles, name) in pred_kcs.items():
        cpds.append((k, smiles, "known", name))
    for k, (smiles, name) in helpers.items():
        cpds.append((k, smiles, "helper", name))
    for k, (smiles, name) in checkpoints.items():
        cpds.append((k, smiles, "checkpoint", name))

    cpds = pl.DataFrame(
        cpds,
        schema=compounds_schema,
        orient="row",
    )

    cpds.write_parquet(
        Path(exp_name) / "compounds.parquet"
    )

    rxns = pl.DataFrame(
        rxn_data,
        schema=expansion_reactions_schema,
        orient="row",
    )

    rxns = rxns.sort(
        "size", descending=True
    ).unique(
        subset=['smarts'], keep="first"
    )
    
    rxns.write_parquet(
        Path(exp_name) / "am_rxns.parquet",
    )
    
    gens = pl.DataFrame(
        {
            "half_expansion": ['forward', 'retro'],
            "generation": [fwd_gens, rev_gens],
        },
        schema=gens_schema,
    )
    gens.write_csv(
        Path(exp_name) / "generations.csv",
        separator=',',
    )

if __name__ == "__main__":
    main()