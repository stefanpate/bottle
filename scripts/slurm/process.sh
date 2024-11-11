#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32GB
#SBATCH -t 24:00:00
#SBATCH --job-name="process"
#SBATCH --output=/home/spn1560/bottle/logs/process/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/process/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
forward=2_steps_ccm_v0_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False
reverse=2_steps_bottle_targets_24_to_None_rules_JN3604IMT_rules_co_metacyc_coreactants_sampled_False_pruned_False

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-OcfYl-TY-py3.12/bin/activate
python -u $scripts_dir/process_expansion.py -f $forward -r $reverse