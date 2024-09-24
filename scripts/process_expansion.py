import pandas as pd
from minedatabase.pickaxe import Pickaxe
from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing import set_start_method

from src.config import filepaths
from src.post_processing import Path, PredictedReaction, KnownReaction, Enzyme, DatabaseEntry, get_path_id
from src.rcmcs import extract_operator_patts, calc_lhs_rcmcs
from src.operator_mapping import expand_paired_cofactors, expand_unpaired_cofactors, standardize_template_map
from src.pickaxe_processing import find_paths, prune_pickaxe
from src.utils import load_json, save_json
from src.chem_draw import draw_rxn_svg
from src.thermo.batch_add_eq_compounds import add_compounds_to_eQ

# Push these into a module?
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
import sqlalchemy
from src.thermo.pickaxe_thermodynamics import PickaxeThermodynamics

if __name__ == '__main__':
    set_start_method("spawn")
    parser = ArgumentParser()
    parser.add_argument("fn", help='Expansion filename w/ extension .pk', type=str)
    parser.add_argument("generations", help="Number of generations run in this expansion", type=int)
    args = parser.parse_args()

    pk_fn = args.fn
    generations = args.generations

    # Set params
    raw_dir = filepaths['expansions']
    pk_path = raw_dir / pk_fn
    k_tautomers = 10 # How many top scoring tautomers to generate for operator mapping
    pre_standardized = False # Predicted reactions assumed pre-standardized

    # Load stored paths
    path_filepath = '../artifacts/processed_expansions/found_paths.json'
    predicted_reactions_filepath = "../artifacts/processed_expansions/predicted_reactions.json"
    known_reactions_filepath = "../artifacts/processed_expansions/known_reactions.json"
    stored_paths = load_json(path_filepath)
    stored_predicted_reactions = load_json(predicted_reactions_filepath)
    stored_known_reactions = load_json(known_reactions_filepath)

    # Read in rules
    rules_dir = "../data/rules"
    rules_fns = ["minimal1224_all_uniprot.tsv", "JN3604IMT_rules.tsv"]
    read_pd = lambda fn : pd.read_csv(f"{rules_dir}/{fn}", sep='\t').set_index("Name").drop(columns='Comments')
    min_rules, imt_rules = [read_pd(fn) for fn in rules_fns]

    # Read in known reactions
    known_reaction_bank = load_json("../data/sprhea/sprhea_240310_v3_mapped_no_subunits.json")
    imt2krs = defaultdict(list)
    min2krs = defaultdict(list)
    for k, v in known_reaction_bank.items():
        min2krs[v['min_rule']].append(k)
        if v['imt_rules']:
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

    # Find paths
    print("Finding paths")
    paths, starters, targets = find_paths(pk, generations)

    # Prune pickaxe output
    pk = prune_pickaxe(pk, paths)
    print(f"Pruned pk object to {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")

    # Add compounds to equilibrator
    add_compounds_to_eQ(pk)

    # Create PredictedReaction objects
    new_predicted_reactions = {}
    for sid, tid in paths.keys():
        for path in paths[(sid, tid)]:        
            for rid in path:
                if rid not in stored_predicted_reactions:
                    pr = PredictedReaction.from_pickaxe(pk, rid)
                    new_predicted_reactions[rid] = pr

    # Add KnownReactions to PredictedReactions
    new_known_reactions = {}
    bad_ops = []
    for id, pr in new_predicted_reactions.items():
        analogues = {}
        srt_imt = sorted(pr.operators, key= lambda x : imt2ct.get(x, 0), reverse=True) # If multiple imt operators, start w/ most common
        for imt in srt_imt:
            min = imt.split('_')[0] # Minify imt operator
            if min2ct.get(min, 0) > 0: # Check if minified imt operator maps known reactions
                
                did_map, aligned_smarts, reaction_center = standardize_template_map(
                    rxn=pr.smarts,
                    rule_row=min_rules.loc[min],
                    smi2paired_cof=smi2paired_cof,
                    smi2unpaired_cof=smi2unpaired_cof,
                    return_rc=True,
                    pre_standardized=pre_standardized,
                )

                if did_map: # Minified op recapitulated imt op
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
                if krid in stored_known_reactions: # Load from stored known reactions
                    kr = KnownReaction.from_dict(stored_known_reactions[krid])
                else: # Create new known reaction from bank
                    bank_kr = known_reaction_bank[krid]
                    
                    # Combine all known reaction operators
                    combined_ops = []
                    if bank_kr['min_rule']:
                        combined_ops.append(bank_kr['min_rule'])
                    if bank_kr['imt_rules']:
                        combined_ops += bank_kr['imt_rules']
                    
                    # Create known reaction object
                    kr = KnownReaction(
                        id=krid,
                        smarts=bank_kr['smarts'],
                        operators=combined_ops,
                        enzymes=[Enzyme.from_dict(e) for e in bank_kr['enzymes']],
                        db_entries=[DatabaseEntry.from_dict({'name': 'rhea', 'id': rhea}) for rhea in bank_kr['rhea_ids']],
                        reaction_center=bank_kr['reaction_center'],
                    )
                    new_known_reactions[krid] = kr # Store in dict of new krs
                    
                # RCMCS
                kr_rcts_rc = [
                    kr.smarts.split('>>')[0].split('.'), # Reactants
                    kr.reaction_center, # Reaction center
                ]
                rcmcs = calc_lhs_rcmcs(pr_rcts_rc, kr_rcts_rc, patts=lhs_patts, norm='max')
                pr.rcmcs[krid] = rcmcs

                analogues[krid] = kr # Append predicted reaction analogues

            pr.analogues = analogues # Add analogues to predicted reaction

    # Connect to compound cache
    with open("../artifacts/eq_uris.uri", "r") as f:
        URI_EQ = f.read().strip("\n")
    
    lcp = LocalCompoundCache()
    lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
    
    # Create pk thermo and eQ objects
    print(f"Getting Thermo Values for {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
    PT = PickaxeThermodynamics(lc=lcp)
    PT.generate_eQ_compound_dict_from_pickaxe(pk=pk)
    PT.generate_eQ_reaction_dict_from_pickaxe(pk=pk)
    
    # Create Path objects
    new_paths = {}
    for sid, tid in paths.keys():
        for path in paths[(sid, tid)]:
            pid = get_path_id(path)

            prs = []
            for rid in path:
                if rid in new_predicted_reactions:
                    prs.append(new_predicted_reactions[rid])
                else:
                    prs.append(PredictedReaction.from_dict(stored_predicted_reactions[rid], stored_known_reactions))

            if pid not in stored_paths:                
                # Calculate path mdf
                mdf_res = PT.calculate_pathway_mdf(reaction_id_list=[pr.id for pr in prs])
                mdf = mdf_res.mdf_value
                dG_opt = {k : v.magnitude for k,v in mdf_res.reaction_energies.items()}
                dG_err = {k : v.magnitude for k,v in mdf_res.uncertainties.items()}
                
                if mdf is None:
                    print(f"Failed mdf for path: {pid}")

                # Add new path
                new_paths[pid] = Path(
                    id=pid,
                    starter=starters[sid],
                    target=targets[tid],
                    reactions=prs,
                    mdf=mdf,
                    dG_opt=dG_opt,
                    dG_err=dG_err,
                    sid=sid,
                    tid=tid,
                )

    # Generate rxn svgs
    for prid, pr in new_predicted_reactions.items():
        pr.image = draw_rxn_svg(pr.smarts, pr.id)

    for krid, kr in new_known_reactions.items():
        kr.image = draw_rxn_svg(kr.smarts, kr.id)

    # Add new to old
    new = [new_known_reactions, new_predicted_reactions, new_paths]
    old = [stored_known_reactions, stored_predicted_reactions, stored_paths]

    for n, o in zip(new, old):
        for id, entry in n.items():
            o[id] = entry.to_dict()

    # Save
    save_json(stored_paths, path_filepath)
    save_json(stored_predicted_reactions, predicted_reactions_filepath)
    save_json(stored_known_reactions, known_reactions_filepath)