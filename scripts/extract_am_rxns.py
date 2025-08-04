import hydra
from omegaconf import DictConfig
from minedatabase.pickaxe import Pickaxe
from pathlib import Path
import polars as pl
from src.schemas import expansion_reactions_schema

@hydra.main(version_base=None, config_path="../conf", config_name="extract_am_rxns")
def main(cfg: DictConfig) -> None:

    if cfg.pk_fn is not None:
        pk_fns = [Path(cfg.pk_fn)]
    else:
        pk_fns = list(Path(cfg.filepaths.raw_expansions).glob("*.pk"))

    for pk_fn in pk_fns:
        
        if (Path.cwd() / f"{pk_fn.stem}").exists():
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

        rxn_data = [
            (">>".join(rxn['am_rxn'].split('>>')[::-1]), list(rxn['Operators'])) 
            for rxn in pk.reactions.values() 
            if 'am_rxn' in rxn
        ]

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