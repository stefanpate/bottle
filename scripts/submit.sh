#!/bin/bash
#SBATCH -A p30041
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 04:00:00
#SBATCH --job-name="2genakg2hopa"
#SBATCH -o ../logs/outlog
#SBATCH -e ../logs/errlog
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u run_pickaxe.py
