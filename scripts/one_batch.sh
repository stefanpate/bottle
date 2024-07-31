#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=20GB
#SBATCH -t 10:00:00
#SBATCH --job-name="sphropmap"
#SBATCH -o ../logs/outlog
#SBATCH -e ../logs/errlog
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
ulimit -c 0
module load python/anaconda3.6
source activate mine
python -u map_operators.py /home/spn1560/bottle/data/rules/minimal1224_all_uniprot.tsv /home/spn1560/bottle/data/sprhea/sprhea_v3_part_19.json /home/spn1560/bottle/artifacts/operator_mapping_sprhea_v3_min_ops_part_19.tsv