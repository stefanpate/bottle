from src.data import Path, PredictedReaction
from src.pickaxe_processing import find_paths, prune_pickaxe
from src.utils import load_json
from src.data import *
from minedatabase.pickaxe import Pickaxe

# Set params
pk_fn = "alpha_ketoglutarate_to_hopa_gen_2_tan_sample_0_n_samples_1000.pk"
generations = 2
path_filepath = '../artifacts/found_paths.json'
# predicted_filepath = "../artifacts/predicted_reactions.json"
# known_filepath = "../artifacts/known_reactions.json"
raw_dir = "/home/stef/bottle/data/raw_expansions"
pk_path = f"{raw_dir}/{pk_fn}"
pruned_path = f"/home/stef/bottle/data/pruned_expansions/{pk_fn}"

# Load raw expansion object
pk = Pickaxe()
pk.load_pickled_pickaxe(pk_path)

# Get pathways
print("Finding paths")
paths, starters, targets = find_paths(pk, generations)

# Prune
pk = prune_pickaxe(pk, paths)
print(f"Saving pruned pk object w/ {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
pk.pickle_pickaxe(pruned_path)

# Create post processing objects
tmp = {}
next_path_id = 0 # TODO: CHANGE
# next_path_id = max(stored_paths.keys()) + 1
for sid, tid in paths.keys():
    for path in paths[(sid, tid)]:
        prs = [PredictedReaction.from_pickaxe(pk, rid) for rid in path]
        tmp[next_path_id] = Path(
            id=next_path_id,
            starter=starters[sid],
            target=targets[tid],
            reactions=prs,
            _sid=sid,
            _tid=tid,
        )
        next_path_id += 1
