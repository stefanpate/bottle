#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10GB
#SBATCH -t 12:00:00
#SBATCH --job-name=sprhmop19
#SBATCH --output=../logs/out/sprhmop19
#SBATCH --error=../logs/error/sprhmop19
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u map_operators.py ../data/rules/minimal1224_all_uniprot.tsv ../data/sprhea/sprhea_v3_part_19.json ../artifacts/operator_mapping_sprhea_v3_min_ops_part_19.tsv