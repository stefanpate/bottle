#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16GB
#SBATCH -t 1:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/home/spn1560/bottle/logs/expand/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/expand/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
starters=pivalic_acid
targets=bottle_targets_24
generations=1
processes=1

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-OcfYl-TY-py3.12/bin/activate
python $scripts_dir/run_pickaxe.py $starters $generations -t $targets -p $processes --prune-to-targets