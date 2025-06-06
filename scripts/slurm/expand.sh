#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 4:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/home/spn1560/bottle/logs/expand/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/expand/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
starters=bottle_targets_24
targets=null
generations=4
processes=50
a_plus_b=True
rules=mechinferred_dt_98_rules_w_coreactants
prune=False

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python3 $scripts_dir/expand.py starters=$starters targets=$targets generations=$generations rules=$rules processes=$processes a_plus_b=$a_plus_b prune_to_targets=$prune