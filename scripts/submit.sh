#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 24:00:00
#SBATCH --job-name="val_kr_uniprot"
#SBATCH -o ../logs/outlog_validation_jni_uniprot
#SBATCH -e ../logs/errlog_validation_jni_uniprot
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u validate_known_rxn_uniprot_mcs.py
