#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48GB
#SBATCH -t 16:00:00
#SBATCH --job-name="lin_path"
#SBATCH --output=/home/spn1560/bottle/logs/pathfinding/out/linear_%A
#SBATCH --error=/home/spn1560/bottle/logs/pathfinding/error/linear_%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu

# Args
script=/home/spn1560/bottle/scripts/linear_pathfinding.py
casp_study=260327_target_3hpa

forward_expansions=(
    2026-03-29/09-35-23 # abf mechinferred_dt_956 2 step
    2026-03-29/09-31-07 # ccm_aa mechinferred_dt_956 2 step
)
retro_expansions=(
    2026-03-30/03-49-57 # 3hpa mechinferred_dt_956 3 step
    2026-03-30/03-49-57 # 3hpa mechinferred_dt_956 3 step
)

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script casp_study=$casp_study forward_expansion=${forward_expansions[0]} retro_expansion=${retro_expansions[0]}
python $script casp_study=$casp_study forward_expansion=${forward_expansions[1]} retro_expansion=${retro_expansions[1]}