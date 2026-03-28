#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=80GB
#SBATCH -t 12:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/home/spn1560/bottle/logs/expand/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/expand/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
starters=3hpa
targets=null
generations=3
processes=50
a_plus_b=False
rules=mechinferred_dt_112_rules_w_coreactants
prune=False
retro=True

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python3 $scripts_dir/expand.py starters=$starters targets=$targets generations=$generations rules=$rules processes=$processes a_plus_b=$a_plus_b prune_to_targets=$prune retro=$retro