from src.data import Path, PredictedReaction, KnownReaction, Enzyme, DatabaseEntry
from src.rcmcs import extract_operator_patts, calc_lhs_rcmcs
from src.operator_mapping import expand_paired_cofactors, expand_unpaired_cofactors, template_map
from src.pickaxe_processing import find_paths, prune_pickaxe
from src.utils import load_json
from minedatabase.pickaxe import Pickaxe
import pandas as pd
from collections import defaultdict

# Set params
pk_fn = "alpha_ketoglutarate_to_hopa_gen_2_tan_sample_0_n_samples_1000.pk"
generations = 2
raw_dir = "/home/stef/bottle/data/raw_expansions"
pk_path = f"{raw_dir}/{pk_fn}"
pruned_path = f"/home/stef/bottle/data/pruned_expansions/{pk_fn}"
k_tautomers = 10

# Load found paths
path_filepath = '../artifacts/found_paths.json'
# TODO: Load
# TODO: next_path_id = max(stored_paths.keys()) + 1

# Read in rules
rules_dir = "/home/stef/bottle/data/rules"
rules_fns = ["minimal1224_all_uniprot.tsv", "JN3604IMT_rules.tsv"]
read_pd = lambda fn : pd.read_csv(f"{rules_dir}/{fn}", sep='\t').set_index("Name").drop(columns='Comments')
min_rules, imt_rules = [read_pd(fn) for fn in rules_fns]

# Read in known reactions
kr_db = load_json("/home/stef/bottle/data/sprhea/sprhea_240310_v3_min_mapped.json")
imt2krs = defaultdict(list)
min2krs = defaultdict(list)
for k, v in kr_db.items():
    min2krs[v['min_rule']].append(k)
    for imt in v['imt_rules']:
        imt2krs[imt].append(k)

min2ct = {k : len(v) for k, v in min2krs.items()}
imt2ct = {k : len(v) for k, v in imt2krs.items()}

# Read in cofactor lookup tables
paired_ref = pd.read_csv('../data/cofactors/paired_cofactors_reference.tsv', sep='\t')
unpaired_ref = pd.read_csv('../data/cofactors/unpaired_cofactors_reference.tsv', sep='\t')
smi2paired_cof = expand_paired_cofactors(paired_ref, k=k_tautomers)
smi2unpaired_cof = expand_unpaired_cofactors(unpaired_ref, k=k_tautomers)

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

# Create PredictedReaction objects
predicted_reactions = {}
for sid, tid in paths.keys():
    for path in paths[(sid, tid)]:        
        for rid in path:
            pr = PredictedReaction.from_pickaxe(pk, rid)
            predicted_reactions[rid] = pr

# Add KnownReactions to PredictedReactions
bad_ops = []
for id, pr in predicted_reactions.items():
    krs = []
    srt_imt = sorted(pr['operators'], key= lambda x : imt2ct.get(x, 0), reverse=True)
    for imt in srt_imt:
        min = imt.split('_')[0] # Minify IMT operator
        if min2ct.get(min, 0) > 0:
            
            did_map, aligned_smarts, reaction_center = template_map(
                rxn=pr.smarts,
                rule_row=min_rules.loc[min],
                smi2paired_cof=smi2paired_cof,
                smi2unpaired_cof=smi2unpaired_cof,
                return_rc=True
            )

            if did_map:
                break
            else:
                bad_ops.append((imt, id))
        
    if did_map:
        pr.smarts = aligned_smarts
        pr.reaction_center = reaction_center
        lhs_patts = extract_operator_patts(min_rules.loc[min, 'SMARTS'], side=0)
        pr_rcts = pr.smarts.split(">>")[0].split('.')
        pr_rcts_rc = [pr_rcts, pr.reaction_center]

        for krid in min2krs[min]:
            kr = kr_db[krid]
            kr_rcts = kr['smarts'].split('>>')[0].split('.')
            kr_rc = kr['reaction_center']
            kr_rcts_rc = [kr_rcts, kr_rc]
            
            rcmcs = calc_lhs_rcmcs(pr_rcts_rc, kr_rcts_rc, patts=lhs_patts, norm='max')

            pr.rcmcs[krid] = rcmcs

            krs.append(
                KnownReaction(
                    id=krid,
                    smarts=kr['smarts'],
                    operators= [kr['min_rule']] + kr['imt_rules'],
                    enzymes=[Enzyme.from_dict(e) for e in kr['enzymes']],
                    db_entries=[DatabaseEntry.from_dict('rhea', rhea) for rhea in kr['rhea_ids']]
                )
            )

        pr.analogues = krs

# Create Path objects
new_paths = {}
next_path_id = 0 # TODO: DELETE
for sid, tid in paths.keys():
    for path in paths[(sid, tid)]:        
        prs = [predicted_reactions[rid] for rid in path]

        new_paths[next_path_id] = Path(
            id=next_path_id,
            starter=starters[sid],
            target=targets[tid],
            reactions=prs,
            _sid=sid,
            _tid=tid,
        )
        next_path_id += 1