#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 3:00:00
#SBATCH --job-name="mech_proc"
#SBATCH --output=/home/spn1560/bottle/logs/mech_proc/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/mech_proc/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
scripts_dir=/home/spn1560/bottle/scripts
processes=50 # Make sure this matches #SBATCH -n
casp_study=bottle25

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python -u $scripts_dir/mechanism_processing.py processes=$processes casp_study=$casp_study