#!/bin/bash

# Define the directory containing the files
directory="/home/stef/quest_data/bottle/results/raw_expansions"

# Filenames
filenames=(
    2mg_to_mvacid_gen_3_tan_sample_1_n_samples_1000.pk
    amino_acids_to_bottle_targets_24_gen_3_tan_sample_1_n_samples_1000_rules_JN3604IMT_rules.pk
    amino_acids_to_bottle_targets_24_gen_3_tan_sample_1_n_samples_1000_rules_JN3604IMT_rules_carbonyl_free.pk
)
generations=(3 3 3)

# Get the length of the arrays
length=${#filenames[@]}

# Loop through the indices of the arrays
for ((i=0; i<length; i++)); do
  # Construct the full path for the current filename
  full_path="${directory}/${filenames[i]}"
  
  message="Running script with: Filename='${full_path}', # generations='${generations[i]}'"
  echo "$message"
  
  python process_expansion.py "$full_path" "${generations[i]}"
done
