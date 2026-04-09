#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=32GB
#SBATCH -t 16:00:00
#SBATCH --job-name="path_proc"
#SBATCH --output=/home/spn1560/bottle/logs/path_proc/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/path_proc/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
processes=10 # Make sure this matches #SBATCH -n
casp_study=260327_target_3hpa

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python -u $scripts_dir/analyze_structures.py processes=$processes casp_study=$casp_study
python -u $scripts_dir/analyze_thermo.py casp_study=$casp_study
python -u $scripts_dir/draw_reactions.py casp_study=$casp_study