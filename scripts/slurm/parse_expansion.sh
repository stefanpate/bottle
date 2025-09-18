#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH -t 0:30:00
#SBATCH --job-name="parse_exp"
#SBATCH --output=/home/spn1560/bottle/logs/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/parse_expansion.py
fwd_exp=3_steps_bottle_targets_25_to_None_rules_mechinferred_dt_04_rules_w_coreactants_co_metacyc_coreactants_sampled_False_pruned_False_aplusb_False.pk
rev_exp=3_steps_bottle_targets_25_to_None_rules_mechinferred_dt_04_rules_w_coreactants_co_metacyc_coreactants_sampled_False_pruned_False_aplusb_False.pk

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script rev_exp=$rev_exp