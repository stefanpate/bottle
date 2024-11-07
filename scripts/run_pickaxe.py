import os
CWD = os.path.abspath("")
os.chdir(CWD)

from src.config import filepaths
from minedatabase.pickaxe import Pickaxe
from minedatabase.rules import metacyc_intermediate
from minedatabase.filters import (
    SimilarityFilter,
    SimilaritySamplingFilter,
)

# Directories and files
st_dir = filepaths['starters_targets']
input_cpds = "alpha_ketoglutarate"
target_cpds = "hopa"
input_cpds_fn = st_dir / f"{input_cpds}.csv"
target_cpds_fn = st_dir / f"{target_cpds}.csv"
rule_list = filepaths['rules'] / 'JN3604IMT_rules_carbonyl_free.tsv'
coreactant_list = filepaths['rules'] /  "metacyc_coreactants_carbonyl_free.tsv"

rules_name = rule_list.stem

# Pickaxe settings
processes = 1 # 50
generations = 1

tani_filter = False
tani_threshold = [0, 0, 0.3, 0.3, 0.3]
increasing_tani = False

tani_sample = True # True
sample_size = 1000 # 1000
weight = None # 5


# save_to = f"/projects/b1039/spn1560/bottle/data/raw_expansions/{input_cpds}_to_{target_cpds}_gen_{generations}_tan_sample_{int(tani_sample)}_n_samples_{sample_size}_rules_{rules_name}.pk"
save_to = "test_expansion.pk"

# _, coreactant_list, rule_name = metacyc_intermediate(
#     fraction_coverage=1
# )

pk = Pickaxe(
        coreactant_list=coreactant_list,
        rule_list=rule_list,
        errors=True,
        quiet=True,
        filter_after_final_gen=True,
    )

pk.load_compound_set(compound_file=input_cpds_fn)

pk.load_targets(target_cpds_fn)


# Apply filters
if tani_filter:
    taniFilter = SimilarityFilter(
        crit_similarity=tani_threshold, increasing_similarity=increasing_tani
    )
    pk.filters.append(taniFilter)

if tani_sample:
    taniSampleFilter = SimilaritySamplingFilter(
        sample_size=sample_size, weight=weight
    )
    pk.filters.append(taniSampleFilter)
    pass

# Transform compounds (the main step)
pk.transform_all(processes, generations)
pk.prune_network_to_targets()

pk.pickle_pickaxe(save_to) # Save results