#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 48:00:00
#SBATCH --job-name="ccm2methylenes"
#SBATCH -o outlog
#SBATCH -e errlog
ulimit -c 0
module load python/anaconda3.6
source activate mine
python run_pickaxe.py
