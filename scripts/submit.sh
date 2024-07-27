#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=16GB
#SBATCH -t 08:00:00
#SBATCH --job-name="rhea_smi_name"
#SBATCH -o ../logs/outlog
#SBATCH -e ../logs/errlog
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u download_rhea.py
