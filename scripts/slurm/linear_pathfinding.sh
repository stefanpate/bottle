#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=0
#SBATCH -t 24:00:00
#SBATCH --job-name="lin_path"
#SBATCH --output=/home/spn1560/bottle/logs/pathfinding/out/linear_%A_%a
#SBATCH --error=/home/spn1560/bottle/logs/pathfinding/error/linear_%A_%a
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.pate@northwestern.edu
#SBATCH --array=0

# Args
script=/home/spn1560/bottle/scripts/linear_pathfinding.py
casp_study=260410_target_3hpa
processes=10 # MAKE SURE THIS MATCHES -n above

forward_expansions=(
    2026-04-16/11-01-33 # pivalic acid mechinferred_dt_956 3 step
)
retro_expansions=(
    2026-03-30/03-49-57 # 3hpa mechinferred_dt_956 3 step
)

# Commands
ulimit -c 0
module purge
source /home/spn1560/.cache/pypoetry/virtualenvs/bottle-jRVXeMfS-py3.12/bin/activate
python $script casp_study=$casp_study processes=$processes forward_expansion=${forward_expansions[$SLURM_ARRAY_TASK_ID]} retro_expansion=${retro_expansions[$SLURM_ARRAY_TASK_ID]}