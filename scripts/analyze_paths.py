from src.rxn_ctr_mcs import *
from src.utils import load_json, rm_atom_map_num
from src.post_processing import *
from rdkit.Chem import AllChem
import pickle
from tqdm import tqdm

starters = 'ccm_v0'
targets = 'mvacid'
generations = 4

expansion_dir = '../data/processed_expansions/'
thermo_dir = '../data/thermo/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000" # Expansion file name

# Load processed expansion
with open(expansion_dir + fn + '.pkl', 'rb') as f:
    pe = pickle.load(f)

pr_am_errors = [] # Track predicted rxn am errors
kr_am_errors = [] # Track known rxn am errors
alignment_issues = [] # Track substrate alignment issues
kekulize_issues = []
norm = 'max atoms' # Normalize MCS atom count by larger molecule

# Populate pred_rxns, known rxn prc-mcs slot

pbar = tqdm(pe.predicted_reactions.items())
for prid, pr in pbar:
    rxn_sma1 = pr.smarts

    # Skip pred reactions that trigger RXNMapper atom mapping errors
    try:
        am_rxn_sma1 = atom_map(rxn_sma1)
    except:
        pr_am_errors.append(prid)
        continue

    a = 0 # Number known rxns analyzed
    for n_kr, kr in enumerate(pr.analogues):
        rxn_sma2 = kr.smarts

        # Catch stoichiometry mismatches stemming from pickaxe, early post-processing
        if tuple([len(elt.split('.')) for elt in rxn_sma2.split('>>')]) != tuple([len(elt.split('.')) for elt in rxn_sma1.split('>>')]):
            print(prid, kr.id, 'stoich_error')
            continue

        # Skip pred reactions that trigger RXNMapper atom mapping errors
        try:
            am_rxn_sma2 = atom_map(rxn_sma2)
        except:
            kr_am_errors.append(kr.id)
            continue

        # Construct reaction objects
        rxns = []
        for elt in [am_rxn_sma1, am_rxn_sma2]:
            temp = AllChem.ReactionFromSmarts(elt, useSmiles=True)
            temp.Initialize()
            rxns.append(temp)

        rc_atoms = [elt.GetReactingAtoms() for elt in rxns] # Get reaction center atom idxs

        # Construct rxn ctr mol objs
        try: # REMOVE after addressing KekulizationException in get_sub_mol
            rcs = []
            for i, t_rxn in enumerate(rxns):
                temp = []
                for j, t_mol in enumerate(t_rxn.GetReactants()):
                    temp.append(get_sub_mol(t_mol, rc_atoms[i][j]))
                rcs.append(temp)
        except:
            kekulize_issues.append((prid, kr.id))
            continue

        # Align substrates of the 2 reactions
        rc_idxs = [] # Each element: (idx for rxn 1, idx for rxn 2)
        remaining = [[i for i in range(len(elt))] for elt in rcs]
        while (len(remaining[0]) > 0) & (len(remaining[1]) > 0):
            idx_pair = align_substrates(rcs, remaining)

            if idx_pair is None:
                break
            else:
                rc_idxs.append(idx_pair)
                remaining[0].remove(idx_pair[0])
                remaining[1].remove(idx_pair[1])

        # Skip if you haven't aligned every reactant pred to known
        if len(rc_idxs) < len(rxn_sma1.split('>>')[0].split('.')):
            alignment_issues.append((prid, kr.id))
            continue

        # For reaction 2 (known reaction) Re-order rcs, rc_atoms,
        # internal order of reactants in the reaction object in rxns
        # and the smarts stored in the known_reactions attribute of the
        # associated predicted reaction

        # Sort reaction 2 rc_idxs by reaction 1 rc_idxs
        rxn_1_rc_idxs, rxn_2_rc_idxs = list(zip(*rc_idxs))
        if rxn_1_rc_idxs != rxn_2_rc_idxs:
            rxn_2_rc_idxs, rxn_1_rc_idxs = sort_x_by_y(rxn_2_rc_idxs, rxn_1_rc_idxs)

            # Re-order atom-mapped smarts string, and then update known_rxns entry
            # with de-atom-mapped version of this string because atom mapper changes
            # reactant order and its this order that rcs, rcatoms, rc_idxs all come from
            am_ro_sma2 = am_rxn_sma2.split('>>')[0].split('.') # Get list of reactant strings
            am_ro_sma2 = '.'.join([am_ro_sma2[elt] for elt in rxn_2_rc_idxs]) # Re-join in new order
            am_rxn_sma2 = am_ro_sma2 + '>>' + am_rxn_sma2.split('>>')[1] # Join with products side

            # Re-construct reaction object from re-ordered, am smarts
            temp = AllChem.ReactionFromSmarts(am_rxn_sma2, useSmiles=True)
            temp.Initialize()
            rxns[1] = temp

            rc_atoms[1] = rxns[1].GetReactingAtoms() # Update rc_atoms
            rcs[1] = [get_sub_mol(elt, rc_atoms[1][i]) for i, elt in enumerate(rxns[1].GetReactants())] # Update rc mol obj
        
        kr.smarts = rm_atom_map_num(am_rxn_sma2) # Update known_reaction entry w/ de-am smarts for consistent ordering in vis
        rxns = align_atom_map_nums(rxns, rcs, rc_atoms)

        # Compute MCS seeded by reaction center
        prc_mcs = get_prc_mcs(rxns, rcs, rc_atoms, norm=norm) 
        pr._mcs_analogues[n_kr][0] = prc_mcs # Update pred_rxns
        
        a += 1 # Count known rxn analyzed
        pr.smarts = rm_atom_map_num(am_rxn_sma1) # Update pred_rxn smarts w/ de-am smarts for consistent ordering in vis

    pbar.set_description(f"Predicted rxn {prid[:5]}. {a / (n_kr + 1):.2f} of {n_kr + 1} analogues successfully processed")

# Thermo

# starters = 'succinate'
# targets = 'mvacid'
# generations = 4
# args = ['-s', f"{starters}", '-t', f"{targets}", '-g', str(generations)]
# command = f"source activate /home/stef/miniconda3/envs/thermo && python /home/stef/pickaxe_thermodynamics/path_mdf.py {' '.join(args)}"
# subprocess.run(command, shell=True)

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