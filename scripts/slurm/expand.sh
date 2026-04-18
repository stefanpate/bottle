#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --mem=0
#SBATCH -t 12:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/home/spn1560/bottle/logs/expand/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/expand/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
starters=pivalic_acid
targets=null
generations=3
processes=25 # Make sure this matches #SBATCH -n
a_plus_b=True
rules=mechinferred_dt_956_rules_w_coreactants
prune=False
retro=False

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python3 $scripts_dir/expand.py starters=$starters targets=$targets generations=$generations rules=$rules processes=$processes a_plus_b=$a_plus_b prune_to_targets=$prune retro=$retro