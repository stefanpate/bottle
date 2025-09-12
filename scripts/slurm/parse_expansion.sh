#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH -t 2:00:00
#SBATCH --job-name="parse_exp"
#SBATCH --output=/home/spn1560/bottle/logs/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/parse_expansion.py
fwd_exp=None
rev_exp=None

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script fwd_exp=$fwd_exp rev_exp=$rev_exp