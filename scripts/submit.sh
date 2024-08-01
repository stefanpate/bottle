#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16GB
#SBATCH -t 5:00:00
#SBATCH --job-name="16spart"
#SBATCH -o ../logs/out/16s
#SBATCH -e ../logs/error/log16s
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u map_operators.py ../data/rules/minimal1224_all_uniprot.tsv ../data/sprhea/sprhea_v3_part_16.json ../artifacts/operator_mapping_sprhea_v3_min_ops_part_16.tsv