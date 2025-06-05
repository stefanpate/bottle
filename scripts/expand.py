import hydra
from omegaconf import DictConfig
from minedatabase.pickaxe import Pickaxe
from minedatabase.filters import SimilaritySamplingFilter
from pathlib import Path

@hydra.main(version_base=None, config_path="../conf", config_name="expand")
def main(cfg: DictConfig) -> None:
    pk = Pickaxe(
        coreactant_list=Path(cfg.filepaths.coreactants) /  f"{cfg.coreactants}.tsv",
        rule_list=Path(cfg.filepaths.rules) / f"{cfg.rules}.tsv",
        errors=True,
        quiet=True,
        filter_after_final_gen=True,
    )

    pk.load_compound_set(compound_file=Path(cfg.filepaths.starters_targets) / f"{cfg.starters}.csv")

    if cfg.a_plus_b:
        # TODO: Re-intro the strict known role matching?
        # known_reactions = load_json(filepaths['data'] / "sprhea" / "sprhea_240310_v3_mapped_no_subunits.json")
        # known_reactions = [{'smarts': v['smarts'], 'rules': v['imt_rules'] if v['imt_rules'] else []} for v in known_reactions.values()]
        # pk.set_starters_as_coreactants(known_reactions=known_reactions)
        pk.set_starters_as_coreactants()

    if cfg.targets:
        pk.load_targets(Path(cfg.filepaths.starters_targets) / f"{cfg.targets}.csv")

    if cfg.do_tani_sample:
        tsf = SimilaritySamplingFilter(sample_size=cfg.sample_size, weight=cfg.weight)
        pk.filters.append(tsf)

    pk.transform_all(cfg.processes, cfg.generations) # Expand

    if cfg.prune_to_targets and cfg.targets:
        pk.prune_network_to_targets()

    fn = (f"{cfg.generations}_steps_{cfg.starters}_to_{cfg.targets}_rules_{cfg.rules}"
          f"_co_{cfg.coreactants}_sampled_{cfg.do_tani_sample}_pruned_{cfg.prune_to_targets}_aplusb_{cfg.a_plus_b}")
    
    pk.pickle_pickaxe(f"{fn}.pk") # Save results

if __name__ == '__main__':
    main()