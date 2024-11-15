import pandas as pd
from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing import set_start_method
from time import perf_counter
from tqdm import tqdm
from src.config import filepaths
from src.post_processing import (
    Expansion,
    Path,
    PredictedReaction,
    KnownReaction,
    Enzyme,
    DatabaseEntry,
    get_path_id,
    realign_pred_rxn_to_rule
) 
from src.rcmcs import extract_operator_patts, calc_lhs_rcmcs
from src.operator_mapping import map_rxn2rule
from src.cheminfo_utils import standardize_smarts_rxn
from src.utils import load_json, save_json
from src.chem_draw import draw_reaction

from src.thermo.batch_add_eq_compounds import add_compounds_to_eQ
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
import sqlalchemy
from src.thermo.pickaxe_thermodynamics import PickaxeThermodynamics

if __name__ == '__main__':
    parser = ArgumentParser(
        description=('Processes forward & reverse expansions singly or in combination.'
                     'Finds paths from starters to targets, stores predicted reactions,'
                     'their known analogues with similarity scores, and thermo calculations.')
    )
    parser.add_argument("stored_dir", help="Subdirectory of results/processed_expansions to read & write processed expansions")
    parser.add_argument("-f", "--forward", default=None, help='Filename of forward expansion (w/o extension)', type=str)
    parser.add_argument("-r", "--reverse", default=None, help='Filename of reverse expansion (w/o extension)', type=str)
    parser.add_argument("--do-thermo", action="store_true", help="Does thermo calculations if provided")
    args = parser.parse_args()

    imt_reverses = load_json(filepaths['rules'] / "jnimt_reverses.json")

    print("Loading expansion")
    pk = Expansion(
        forward=filepaths['raw_expansions'] / f"{args.forward}.pk" if args.forward else args.forward,
        reverse=filepaths['raw_expansions'] / f"{args.reverse}.pk" if args.reverse else args.reverse,
        operator_reverses=imt_reverses,
    )
    print("Searching for paths")
    tic = perf_counter()
    paths = pk.find_paths()
    toc = perf_counter()
    print(f"Path finding took: {toc - tic : .2f} seconds")
    pk.prune(paths)
    print(f"Pruned expansion to {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")

    if args.do_thermo:
        set_start_method("spawn")
    
    # Set params
    pre_standardized = False # Predicted reactions assumed pre-standardized

    # Load stored paths
    stored_fp = filepaths['processed_expansions'] / args.stored_dir 
    if not stored_fp.exists():
        stored_fp.mkdir()

    path_filepath = stored_fp / 'found_paths.json'
    predicted_reactions_filepath = stored_fp / "predicted_reactions.json"
    known_reactions_filepath = stored_fp / "known_reactions.json"
    load_processed = lambda path : load_json(path) if path.exists() else {}
    stored_paths = load_processed(path_filepath)
    stored_predicted_reactions = load_processed(predicted_reactions_filepath)
    stored_known_reactions = load_processed(known_reactions_filepath)

    # Read in rules
    rules_dir = filepaths['rules']
    rules_fns = ["minimal1224_all_uniprot.tsv", "JN3604IMT_rules.tsv"]
    read_pd = lambda fn : pd.read_csv(f"{rules_dir}/{fn}", sep='\t').set_index("Name").drop(columns='Comments')
    min_rules, imt_rules = [read_pd(fn) for fn in rules_fns]

    # Read in known reactions
    # NOTE: Minified IMT operator must match MIN operator of a known reaction to be counted
    known_reaction_bank = load_json(filepaths['data'] / "sprhea/sprhea_240310_v3_mapped_no_subunits.json")
    imt2krs = defaultdict(list)
    for k, v in known_reaction_bank.items():
        if v['imt_rules'] and v['min_rule']:
            for imt in v['imt_rules']:
                if imt.split('_')[0] == v['min_rule']:
                    imt2krs[imt].append(k)

    imt2ct = {k : len(v) for k, v in imt2krs.items()}

    if args.do_thermo:
        print("Adding compounds to equilibrator")
        add_compounds_to_eQ(pk)

    # Create new PredictedReaction objects where don't already have
    print("Creating new predicted reactions")
    new_predicted_reactions = {}
    for sid, tid in paths.keys():
        for path in paths[(sid, tid)]:        
            for rid in path:
                if rid not in stored_predicted_reactions:
                    pr = PredictedReaction.from_pickaxe(pk, rid)
                    new_predicted_reactions[rid] = pr

    # Add KnownReactions to new PredictedReactions
    print("Adding known analogues")
    tic = perf_counter()
    new_known_reactions = {}
    unmapped = dict(new_predicted_reactions.items())
    for id, pr in tqdm(sorted(new_predicted_reactions.items())):
        most_common_map = None
        analogues = {}
        srt_imt = sorted(pr.operators, key= lambda x : imt2ct.get(x, 0), reverse=True) # If multiple imt operators, start w/ most common
        for imt in srt_imt:
            min = imt.split('_')[0] # Minify imt operator to get reaction center by protection-guess-and-check

            # Standardize smarts
            if not pre_standardized:
                try:
                    rxn = standardize_smarts_rxn(pr.smarts, quiet=True)
                except:
                    print(f"Unable to standardize reaction: {pr.smarts}")
                    rxn = pr.smarts

            # TODO things could get tricky here with combo expansions and multiple sets of coreactants...
            matched_idxs = realign_pred_rxn_to_rule(rxn, min_rules.loc[min, "Reactants"], pk.coreactants) # Align reactants to rule
                
            # Map rule to reaction
            if matched_idxs:
                res = map_rxn2rule(rxn, min_rules.loc[min, "SMARTS"], return_rc=True, matched_idxs=matched_idxs)
                did_map, aligned_smarts, reaction_center = res['did_map'], res['aligned_smarts'], res['reaction_center']
            else:
                did_map = False

            if did_map:

                if most_common_map is None: # First operator mapped is the most common
                    pr.smarts = aligned_smarts
                    pr.reaction_center = reaction_center
                    lhs_patts = extract_operator_patts(min_rules.loc[min, 'SMARTS'], side=0)
                    pr_rcts = pr.smarts.split(">>")[0].split('.')
                    pr_rcts_rc = [pr_rcts, pr.reaction_center]
                    unmapped.pop(id, None)
                    most_common_map = imt
                
                # If this is the first op or subsequent but same underlying min rule
                # and therefore reaction center, then you can assign analogues
                if most_common_map.split('_')[0] == imt.split('_')[0]:
                    for krid in imt2krs[imt]:
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
            else:
                print(f"Minified operator failed to recapitulate reaction {imt} {id}")

            pr.analogues = analogues # Add analogues to predicted reaction

    toc = perf_counter()
    print(f"Analogue analysis took: {toc - tic : .2f} seconds")
    print(f"{len(unmapped)} predicted reactions left unmapped")
    
    if args.do_thermo:
        # Connect to compound cache
        with open(filepaths['artifacts'] / "eq_uris.uri", "r") as f:
            URI_EQ = f.read().strip("\n")
        
        lcp = LocalCompoundCache()
        lcp.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
    
        # Create pk thermo and eQ objects
        print(f"Getting Thermo Values for {len(pk.compounds)} compounds and {len(pk.reactions)} reactions")
        PT = PickaxeThermodynamics(lc=lcp)
        PT.generate_eQ_compound_dict_from_pickaxe(pk=pk)
        PT.generate_eQ_reaction_dict_from_pickaxe(pk=pk)
    
    # Create Path objects
    print("Adding new paths (and calculating mdf)")
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

            # If new path or missing thermo, and want to do thermo now, do it
            if (pid not in stored_paths or stored_paths[pid]['mdf'] is None) and args.do_thermo:
                mdf_res = PT.calculate_pathway_mdf(reaction_id_list=[pr.id for pr in prs])
                mdf = mdf_res.mdf_value
                dG_opt = {k : v.magnitude for k,v in mdf_res.reaction_energies.items()}
                dG_err = {k : v.magnitude for k,v in mdf_res.uncertainties.items()}
                
                if mdf is None:
                    print(f"Failed mdf for path: {pid}")
            
            # Otherwise these are the defaults
            else:
                mdf = None
                dG_opt = {}
                dG_err = {}

            # If it was new, create a new path
            if pid not in stored_paths:
                # Add new path
                new_paths[pid] = Path(
                    id=pid,
                    starter=pk.starters[sid],
                    target=pk.targets[tid],
                    reactions=prs,
                    mdf=mdf,
                    dG_opt=dG_opt,
                    dG_err=dG_err,
                    sid=sid,
                    tid=tid,
                )
            
            # If it was missing thermo and old, update stored paths
            elif stored_paths[pid]['mdf'] is None and args.do_thermo:
                stored_paths[pid]['mdf'] = mdf
                stored_paths[pid]['dG_opt'] = dG_opt
                stored_paths[pid]['dG_err'] = dG_err

    # TODO: Align w/ new chem_draw api
    # Generate rxn svgs
    for prid, pr in new_predicted_reactions.items():
        pr.image = draw_reaction(pr.smarts, pr.id, auto_scl=True)

    for krid, kr in new_known_reactions.items():
        kr.image = draw_reaction(kr.smarts, kr.id, auto_scl=True)

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