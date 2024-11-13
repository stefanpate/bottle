#!/bin/bash

stored_dir=241118_bottle_checkin
proc_path=/home/stef/bottle/scripts/process_expansion.py
python $proc_path $stored_dir -f 1_steps_dmb_to_dmhb_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_True --do-thermo
python $proc_path $stored_dir -f 2_steps_pivalic_acid_to_bottle_targets_24_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_True --do-thermo
python $proc_path $stored_dir -f 2_steps_ccm_v0_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False -r 2_steps_bottle_targets_24_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False --do-thermo