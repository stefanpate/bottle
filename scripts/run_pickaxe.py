import os
CWD = os.path.abspath("")
os.chdir(CWD)

from minedatabase.pickaxe import Pickaxe
from minedatabase.utils import get_compound_hash
from minedatabase.rules import metacyc_intermediate
from minedatabase.filters import (
    SimilarityFilter,
    SimilaritySamplingFilter,
)

from rdkit.Chem import CanonSmiles

# Directories and files
st_dir = "../src/starters_targets/"
input_cpds = st_dir + "methylene_molecules.csv"
target_cpds = st_dir + "mvacid.csv"
rule_list = '../src/rules/JN3604IMT_rules.tsv'

# Pickaxe settings
processes = 50
generations = 4

tani_filter = False # True
tani_threshold = [0, 0, 0.3, 0.3, 0.3]
increasing_tani = False

tani_sample = True # False
sample_size = 1000 # 10
weight = None # 5

save_to = f"/projects/b1039/spn1560/bottle/{input_cpds.split('/')[-1]}_to_{target_cpds.split('/')[-1]}_gen_{generations}_tan_sample_{int(tani_sample)}_n_samples_{sample_size}.pk"


_, coreactant_list, rule_name = metacyc_intermediate(
    fraction_coverage=1
    # n_rules=5
)

pk = Pickaxe(
        coreactant_list=coreactant_list,
        rule_list=rule_list,
        errors=True,
        quiet=True,
        filter_after_final_gen=True,
    )

pk.load_compound_set(compound_file=input_cpds)

pk.load_targets(target_cpds)


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
