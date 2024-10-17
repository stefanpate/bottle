#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 36:00:00
#SBATCH --job-name="map"
#SBATCH -o ../logs/out/min_remap_v3
#SBATCH -e ../logs/error/e_tmp_1
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate casp
python -u map_operators.py minimal1224_all_uniprot.tsv sprhea/sprhea_240310_v3.json operator_mapping_sprhea_v3_min_ops_new.tsv