from src.rxn_ctr_mcs import *
from src.utils import load_json, rm_atom_map_num
from src.post_processing import *

from minedatabase.pickaxe import Pickaxe

from rdkit.Chem import AllChem
import pickle

# Load processed expansion
starters = '2mg'
targets = 'mvacid'
generations = 2

expansion_dir = '../data/processed_expansions/'
thermo_dir = '../data/thermo/'
fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name

# Load processed expansion
with open(expansion_dir + fn, 'rb') as f:
    pe = pickle.load(f)

pr_am_errors = [] # Track predicted rxn am errors
kr_am_errors = [] # Track known rxn am errors
alignment_issues = [] # Track substrate alignment issues
kekulize_issues = []
norm = 'max atoms' # Normalize MCS atom count by larger molecule

# Populate pred_rxns, known rxn prc-mcs slot
for prid, pr in pe.predicted_reactions.items():
    rxn_sma1 = pr.smarts

    # Skip pred reactions that trigger RXNMapper atom mapping errors
    try:
        am_rxn_sma1 = atom_map(rxn_sma1)
    except:
        pr_am_errors.append(prid)
        continue

    a = 0 # Number known rxns analyzed
    for z, kr in enumerate(pr.analogues):
        rxn_sma2 = kr.smarts

        # Catch stoichiometry mismatches stemming from pickaxe, early post-processing
        if tuple([len(elt.split('.')) for elt in rxn_sma2.split('>>')]) != tuple([len(elt.split('.')) for elt in rxn_sma1.split('>>')]):
            print(prid, z, 'stoich_error')
            continue

        # Skip pred reactions that trigger RXNMapper atom mapping errors
        try:
            am_rxn_sma2 = atom_map(rxn_sma2)
        except:
            kr_am_errors.append((prid, j))
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
            kekulize_issues.append((prid, z))
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
            alignment_issues.append((prid, z))
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
        
        kr.smarts = rm_atom_map_num(am_rxn_sma2) # Update known_reaction entry w/ de-am smarts
        rxns = align_atom_map_nums(rxns, rcs, rc_atoms)

        # Compute MCS seeded by reaction center
        prc_mcs = get_prc_mcs(rxns, rcs, rc_atoms, norm=norm) 
        pr.prc_mcs = prc_mcs # Update pred_rxns
        
        a += 1 # Count known rxn analyzed
        pr.smarts = rm_atom_map_num(am_rxn_sma1) # Update pred_rxn smarts w/ de-am smarts

    print(prid[:5], ':', a / (z+1), 'of', z+1)

# Thermo

# starters = 'succinate'
# targets = 'mvacid'
# generations = 4
# args = ['-s', f"{starters}", '-t', f"{targets}", '-g', str(generations)]
# command = f"source activate /home/stef/miniconda3/envs/thermo && python /home/stef/pickaxe_thermodynamics/path_mdf.py {' '.join(args)}"
# subprocess.run(command, shell=True)

#TODO: Check how themo json being generated. Need assurance that the order
# of this file lines up with paths in pe object
thermo = load_json(thermo_dir + fn)
for k,v in thermo.items():
    st = tuple(k.split('>'))
    for i, elt in enumerate(thermo[k]):
        paths[st][i].mdf = elt['mdf']
        paths[st][i].dG_opt = elt['dG_opt']
        paths[st][i].dG_err = elt['dG_err']
        paths[st][i].conc_opt = elt['conc_opt']


# Save reactions dict and paths list (ultimately will replace with expansion object)

rxns_fn = 'predicted_reactions_' + fn
paths_fn = 'paths_' + fn
save_dir = '../data/processed_expansions/'
rxns_path = save_dir + rxns_fn
paths_path = save_dir + paths_fn

with open(rxns_path, 'wb') as f:
    pickle.dump(pred_rxns, f)

with open(paths_path, 'wb') as f:
    pickle.dump(paths, f)