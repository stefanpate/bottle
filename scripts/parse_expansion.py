import hydra
from omegaconf import DictConfig
from minedatabase.pickaxe import Pickaxe
from pathlib import Path
import polars as pl
from src.schemas import expansion_reactions_schema
from tqdm import tqdm
from itertools import chain

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

@hydra.main(version_base=None, config_path="../conf", config_name="parse_expansion")
def main(cfg: DictConfig) -> None:

    if cfg.pk_fn is not None:
        pk_fns = [Path(cfg.pk_fn)]
    else:
        pk_fns = list(Path(cfg.filepaths.raw_expansions).glob("*.pk"))

    for pk_fn in pk_fns:
        
        if Path(f"{pk_fn.stem}").exists():
            continue
        else:
            Path(f"{pk_fn.stem}").mkdir()
        
        pk = Pickaxe()
        pk.load_pickled_pickaxe(
            Path(cfg.filepaths.raw_expansions) / pk_fn
        )

        if cfg.mode == "retro":
            kcs = pl.read_parquet(
                Path(cfg.filepaths.known) / "known_compounds.parquet",
            )

            krs = pl.read_parquet(
                Path(cfg.filepaths.known) / "known_reactions.parquet",
            )
            reaction_reverses = dict(zip(krs["id"], krs["reverse"]))

            mapped_rxns = pl.read_parquet(
                Path(cfg.rxn_x_rule_mapping),
            )
            rule2rxn = mapped_rxns.group_by("rule_id").agg(
                pl.col("rxn_id")
            )
            rule2rxn = dict(zip(rule2rxn["rule_id"], rule2rxn["rxn_id"].to_list()))
            rxn2rule = dict(zip(mapped_rxns["rxn_id"], mapped_rxns["rule_id"].to_list()))
            targets = [cpd['SMILES'] for cpd in pk.compounds.values() if cpd['Type'] == 'Starting Compound']

            pred_kcs = []
            kcs_smiles = set(kcs["smiles"].to_list())
            for cid, cpd in pk.compounds.items():
                if cpd['SMILES'] in kcs_smiles and cpd['Type'] != 'Starting Compound':
                    pred_kcs.append(cpd['SMILES'])

            with open(Path(pk_fn.stem) / "pred_kcs.txt", "w") as f:
                for smi in pred_kcs:
                    f.write(smi + "\n")

            with open(Path(pk_fn.stem) / "targets.txt", "w") as f:
                for smi in targets:
                    f.write(smi + "\n")

        rxn_data = []
        for rxn in tqdm(pk.reactions.values(), desc="Processing reactions", total=len(pk.reactions)):
            if 'am_rxn' not in rxn:
                continue
            
            rules = [elt.split('_')[0] for elt in rxn['Operators']]
            
            if cfg.mode == "retro":
                am_smarts = ">>".join(rxn['am_rxn'].split('>>')[::-1])
                rules = reverse_rules(rules, rule2rxn, rxn2rule, reaction_reverses)
            else:
                am_smarts = rxn['am_rxn']

            rules = [f"{cfg.rule_set}:{rule}" for rule in rules]
            rxn_data.append((am_smarts, rules))

        rxns = pl.DataFrame(
            rxn_data,
            schema=expansion_reactions_schema,
            orient="row",
        )

        rxns.write_parquet(
            Path(pk_fn.stem) / "am_rxns.parquet",
        )

if __name__ == "__main__":
    main()