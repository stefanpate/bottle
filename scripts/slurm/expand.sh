#!/bin/bash
#SBATCH -A b1039
#SBATCH -p b1039
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mem=0
#SBATCH -t 4:00:00
#SBATCH --job-name="exp"
#SBATCH --output=/path/to/outlogs/%A
#SBATCH --error=/path/to/error/logs/%A
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<YOUR EMAIL>

# Args
scripts_dir=/path/to/your/scripts
starters=foo
generations=1
targets=bar

# Commands
ulimit -c 0
module purge
source /path/to/your/environment/bin/activate
python $scripts_dir/run_pickaxe.py $starters $generations