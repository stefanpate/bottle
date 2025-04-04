from argparse import ArgumentParser
from src.config import filepaths
from src.utils import load_json
from src.feasfilter import DORAXGBFilter
from minedatabase.pickaxe import Pickaxe
from minedatabase.filters import SimilaritySamplingFilter

def main(args):
    pk = Pickaxe(
        coreactant_list=filepaths['coreactants'] /  f"{args.coreactants}.tsv",
        rule_list=filepaths['rules'] / f"{args.rules}.tsv",
        errors=True,
        quiet=True,
        filter_after_final_gen=True,
    )

    pk.load_compound_set(compound_file=filepaths['starters_targets'] / f"{args.starters}.csv")

    if args.a_plus_b:
        known_reactions = load_json(filepaths['data'] / "sprhea" / "sprhea_240310_v3_mapped_no_subunits.json")
        known_reactions = [{'smarts': v['smarts'], 'rules': v['imt_rules'] if v['imt_rules'] else []} for v in known_reactions.values()]
        pk.set_starters_as_coreactants(known_reactions=known_reactions)

    if args.targets:
        pk.load_targets(filepaths['starters_targets'] / f"{args.targets}.csv")

    if args.tani_sample:
        tsf = SimilaritySamplingFilter(sample_size=args.sample_size, weight=args.weight)
        pk.filters.append(tsf)

    if args.feas_filter:  # Corrected the undefined variable
        feasibility_filter = DORAXGBFilter(threshold=args.feas_threshold, generation_list=None, last_generation_only=False)
        pk.filters.append(feasibility_filter)

    pk.transform_all(args.processes, args.generations) # Expand

    if args.prune_to_targets and args.targets:
        pk.prune_network_to_targets()

    fn = (f"{args.generations}_steps_{args.starters}_to_{args.targets}_rules_{args.rules}"
          f"_co_{args.coreactants}_sampled_{args.tani_sample}_pruned_{args.prune_to_targets}_aplusb_{args.a_plus_b}")
    
    pk.pickle_pickaxe(filepaths['raw_expansions'] / f"{fn}.pk") # Save results

parser = ArgumentParser(description="Run expansion")
parser.add_argument("starters", help="Name of starters file w/o extension")
parser.add_argument("generations", type=int, help="Number of times to expand")
parser.add_argument("-t", "--targets", default=None, help="Name of targets file w/o extension")
parser.add_argument("-r", "--rules", default="JN3604IMT_rules", help="Name of rules file w/o extension")
parser.add_argument("-c", "--coreactants", default="metacyc_coreactants", help="Name of coreactants file w/o extension")
parser.add_argument("-p", "--processes", type=int, default=1, help="Number of cpus")
parser.add_argument("-s", "--sample-size", type=int, default=1000, help="Number of compounds to sample each generation")
parser.add_argument("-w", "--weight", type=int, default=None, help="Sampling distribution parameter")
parser.add_argument("--prune-to-targets", action="store_true", help="Keep only compounds that lead to target")
parser.add_argument("--tani-sample", action="store_true", help="Sample compounds as fcn of tanimoto sim to targets")
parser.add_argument("--a-plus-b", action='store_true', help="Allow use of starters as coreactants in multisubstrate reactions")
parser.add_argument("--feas-filter", action="store_true", help="Apply feasibility filter")  # Added feasibility filter option
parser.add_argument("--feas-threshold", type=float, default=0.5, help="Threshold for feasibility filtering")  # Added feasibility threshold

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)