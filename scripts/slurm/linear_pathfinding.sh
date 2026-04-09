#!/bin/bash
#SBATCH -A p30041
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24GB
#SBATCH -t 4:00:00
#SBATCH --job-name="lin_path"
#SBATCH --output=/home/spn1560/bottle/logs/out/%A
#SBATCH --error=/home/spn1560/bottle/logs/error/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/linear_pathfinding.py
forward_expansion=2026-03-26/17-52-16
retro_expansion=2026-03-26/17-54-12
casp_study=260327_target_3hpa
processes=10

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script forward_expansion=$forward_expansion retro_expansion=$retro_expansion casp_study=$casp_study processes=$processes