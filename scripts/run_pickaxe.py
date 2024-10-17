import os
CWD = os.path.abspath("")
os.chdir(CWD)

from minedatabase.pickaxe import Pickaxe
from minedatabase.rules import metacyc_intermediate
from minedatabase.filters import (
    SimilarityFilter,
    SimilaritySamplingFilter,
)

# Directories and files
st_dir = "../data/starters_targets/"
input_cpds = "amino_acids"
target_cpds = "bottle_targets_24"
input_cpds_fn = st_dir + input_cpds + ".csv"
target_cpds_fn = st_dir + target_cpds + ".csv"
rule_list = '../data/rules/JN3604IMT_rules_carbonyl_free.tsv'
coreactant_list = "../data/rules/metacyc_coreactants_carbonyl_free.tsv"

rules_name = rule_list.split('/')[-1][:-4]

# Pickaxe settings
processes = 50 # 50
generations = 3
processes = 50 # 50
generations = 3

tani_filter = False
tani_threshold = [0, 0, 0.3, 0.3, 0.3]
increasing_tani = False

tani_sample = True # True
sample_size = 1000 # 1000
weight = None # 5


save_to = f"/projects/b1039/spn1560/bottle/data/raw_expansions/{input_cpds}_to_{target_cpds}_gen_{generations}_tan_sample_{int(tani_sample)}_n_samples_{sample_size}_rules_{rules_name}.pk"

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