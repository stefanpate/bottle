from src.rxn_ctr_mcs import *
from src.data import *
import pickle
from tqdm import tqdm
import pandas as pd

starters = 'alpha_ketoglutarate'
targets = 'hopa'
generations = 2
ts = 0

expansion_dir = '../data/processed_expansions/'
thermo_dir = '../data/thermo/'
norm = 'max atoms' # Normalize MCS atom count by larger molecule
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_{ts}_n_samples_1000" # Expansion file name

# Load processed expansion
with open(expansion_dir + fn + '.pkl', 'rb') as f:
    pe = pickle.load(f)

min_rules = pd.read_csv("../data/rules/minimal1224_all_uniprot.tsv", sep='\t')
min_rules.set_index("Name", inplace=True)

# Populate pred_rxns, known rxn prc-mcs slot
pbar = tqdm(pe.predicted_reactions.items())
matched_krs = set()
for prid, pr in pbar:
    potential_perm_idxs = defaultdict(list) # To store possible perm idx, mcs where mutliple between PR & KR are possible
    
    # Get reaction center of LHS of PR
    pr_smarts = pr.smarts
    pr_rcs = []
    for imt_rule in pr.imt_rules:
        min_rule_name = imt_rule.split('_')[0]
        min_rule = min_rules.loc[min_rule_name, 'SMARTS']
        pr_rc, pr_smarts_perm = get_pred_rxn_ctr(pr_smarts, min_rule) # PR RC + permuted LHS of smarts to match op template order

        if pr_rc is None:
            print(pr.id, imt_rule, pr.imt_rules)
            continue

        pr_rcs.append(pr_rc)

        for n_kr, kr in enumerate(pr.analogues):
            kr_smarts = kr.smarts
            for i, kr_rc in enumerate(kr.rcs):
                kr_min_rule_name = kr.min_rules[i]
                
                if kr_min_rule_name != min_rule_name:
                    print(f"Rule mistmatch")
                    continue # PR and KR min rule names did not match. TODO: eliminate these from PE?

                if tuple([len(elt.split('.')) for elt in pr_smarts.split('>>')]) != tuple([len(elt.split('.')) for elt in kr_smarts.split('>>')]):
                    print(f"Stoich error")
                    continue

                # Returns permuted indices of PR to align to KR RCs
                kr_smarts_perm, kr_rc_perm = align_reactants_to_operator(
                    reaction_smarts=kr_smarts,
                    reaction_center=kr_rc,
                    min_rule=min_rule
                    )

                matched_krs.add(kr.id)

                # Compute rc-mcs values
                rxns = [pr_smarts_perm, kr_smarts_perm]
                rc_atoms = [pr_rc, kr_rc_perm]
                rc_mcs = get_rc_mcs(
                    rxns,
                    rc_atoms,
                    min_rule,
                    norm=norm)
                
                potential_perm_idxs[n_kr].append((kr_smarts_perm, rc_mcs)) # TODO: also take rc perm and  min & imt operators to update kr entry in pe

    # Take perm_idx of the best rc_mcs
    pr.smarts = pr_smarts_perm # TODO: need to have multiple PR entries now in case of different perms?
    for k, v in potential_perm_idxs.items():
        print(f"{len(potential_perm_idxs) / len(pr.analogues)} of {len(pr.analogues)} analogues processed for {pr.id[:5]}")
        best_entry = sorted(v, key=lambda x : sum(x[1]) / len(x[1]), reverse=True)[0]
        best_smarts, best_rc_mcs = best_entry
        pr._mcs_analogues[n_kr][0] = best_rc_mcs
        pr._mcs_analogues[n_kr][1].smarts = best_smarts
        # print(pr.id, '\n', pr.smarts, '\n', best_smarts, '\n', best_rc_mcs)

not_matched_krs = set()
for prid, pr in pbar:    
    for n_kr, kr in enumerate(pr.analogues):
        if kr.id not in matched_krs:
            not_matched_krs.add(kr.id)
            
# Thermo
# TODO: bring thermo within this directory
# starters = 'succinate'
# targets = 'mvacid'
# generations = 4
# args = ['-s', f"{starters}", '-t', f"{targets}", '-g', str(generations)]
# command = f"source activate /home/stef/miniconda3/envs/thermo && python /home/stef/pickaxe_thermodynamics/path_mdf.py {' '.join(args)}"
# subprocess.run(command, shell=True)

# Get thermo calcs from pickaxe_thermodynamics
thermo = load_json(thermo_dir + fn + '.json')
for k,v in thermo.items():
    st = tuple(k.split('>'))
    for i, elt in enumerate(thermo[k]):
        if elt['mdf']:
            pe._st2paths[st][i].mdf = elt['mdf']
            pe._st2paths[st][i].dG_opt = elt['dG_opt']
            pe._st2paths[st][i].dG_err = elt['dG_err']
            pe._st2paths[st][i].conc_opt = elt['conc_opt']


# Save processed expansion object
print("Saving processed expansion object")
save_dir = '../data/processed_expansions/'
with open(save_dir + fn + '.pkl', 'wb') as f:
    pickle.dump(pe, f)