import hydra
from omegaconf import DictConfig
from minedatabase.pickaxe import Pickaxe
from pathlib import Path
import polars as pl
from src.schemas import expansion_reactions_schema
from tqdm import tqdm
from itertools import chain
from src.network import de_am
from ergochemics.mapping import rc_to_nest

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
        
        rules = [elt.split('_')[0] for elt in rxn['Operators']]
        
        if mode == "retro":
            am_smarts = ">>".join(rxn['am_rxn'].split('>>')[::-1])
            rules = reverse_rules(rules, rule2rxn, rxn2rule, reaction_reverses)
        else:
            am_smarts = rxn['am_rxn']

        lhs, rhs = de_am(am_smarts)
        smarts = f"{".".join(lhs)}>>{".".join(rhs)}"
        size = max([rule2size.get(rule, 0) for rule in rules])
        rules = [f"{rule_set}:{rule}" for rule in rules]
        rxn_data.append((smarts, am_smarts, rules, size))

    return rxn_data

@hydra.main(version_base=None, config_path="../conf", config_name="parse_expansion")
def main(cfg: DictConfig) -> None:
    
    if cfg.fwd_exp is not None and cfg.rev_exp is not None:
        steps = int(cfg.fwd_exp[0]) + int(cfg.rev_exp[0])
        fwd_rules = cfg.fwd_exp.split('_rules_')[1]
        rev_rules = cfg.rev_exp.split('_rules_')[1]
        exp_name = Path(f"{steps}_steps_combo_{fwd_rules}_and_{rev_rules}_rules")
    elif cfg.fwd_exp is not None:
        exp_name = Path(cfg.fwd_exp).stem
    elif rev_exp is not None:
        exp_name = Path(cfg.rev_exp).stem
  
    if not Path(exp_name).exists():
        Path(exp_name).mkdir(parents=True, exist_ok=True)
    
    rxn_data = []
    helpers = set()
    if cfg.fwd_exp is not None:
        fwd_rule_set = cfg.fwd_exp.split('_rules_')[1]
        fwd_exp = Pickaxe()
        fwd_exp.load_pickled_pickaxe(
            Path(cfg.filepaths.raw_expansions) / Path(cfg.fwd_exp)
        )

        mapped_rxns = pl.read_parquet(
            Path(cfg.filepaths.rxn_x_rule_mapping) / f"mapped_known_reactions_x_{fwd_rule_set}_rules.parquet",
        )

        rule2size = mapped_rxns.group_by("rule_id").agg(
            pl.col("template_aidxs").first().map_elements(lambda x: sum(len(elt) for elt in rc_to_nest(x)[0]), return_dtype=pl.Int32).alias("size")
        )

        rule2size = dict(zip(rule2size["rule_id"], rule2size["size"].to_list()))

        sources = []
        for cid, cpd in fwd_exp.compounds.items():
            if cid[0] == 'X':
                helpers.add(cpd['SMILES'])

            if cpd['Type'] == 'Starting Compound':
                sources.append(cpd['SMILES'])
            
        with open(Path(exp_name) / "sources.txt", "w") as f:
            for smi in sources:
                f.write(smi + "\n")

        targets = [cpd["SMILES"] for cpd in fwd_exp.targets.values()]

        with open(Path(exp_name) / "targets.txt", "w") as f:
            for smi in targets:
                f.write(smi + "\n")

        rxn_data += parse_reactions(
            fwd_exp.reactions,
            rule2rxn={},
            rxn2rule={},
            rule2size=rule2size,
            reaction_reverses={},
            mode="fwd",
            rule_set=fwd_rule_set,
        )

        del fwd_exp
  
    if cfg.rev_exp is not None:
        rev_rule_set = cfg.rev_exp.split('_rules_')[1]
        rev_exp = Pickaxe()
        rev_exp.load_pickled_pickaxe(
            Path(cfg.filepaths.raw_expansions) / Path(cfg.rev_exp)
        )

        kcs = pl.read_parquet(
            Path(cfg.filepaths.known) / "known_compounds.parquet",
        )

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
        targets = [cpd['SMILES'] for cpd in rev_exp.compounds.values() if cpd['Type'] == 'Starting Compound']

        pred_kcs = []
        kcs_smiles = set(kcs["smiles"].to_list())
        for cid, cpd in rev_exp.compounds.items():
            if cpd['SMILES'] in kcs_smiles and cpd['Type'] != 'Starting Compound':
                pred_kcs.append(cpd['SMILES'])

            if cid[0] == 'X':
                helpers.add(cpd['SMILES'])

        with open(Path(exp_name) / "pred_kcs.txt", "w") as f:
            for smi in pred_kcs:
                f.write(smi + "\n")

        with open(Path(exp_name) / "targets.txt", "w") as f:
            for smi in targets:
                f.write(smi + "\n")

        rxn_data += parse_reactions(
            rev_exp.reactions,
            rule2rxn=rule2rxn,
            rxn2rule=rxn2rule,
            rule2size=rule2size,
            reaction_reverses=reaction_reverses,
            mode="retro",
            rule_set=rev_rule_set,
        )

    with open(Path(exp_name) / "helpers.txt", "w") as f:
        for smi in helpers:
            f.write(smi + "\n")

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

if __name__ == "__main__":
    main()