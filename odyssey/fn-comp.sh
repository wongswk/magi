#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 26 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 16000 # Memory request
#SBATCH -t 2-01:20 # (D-HH:MM)
#SBATCH -o /n/scratchlfs/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.out # Standard output
#SBATCH -e /n/scratchlfs/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.err # Standard error


source "$HOME/R/load.sh"
source "$HOME/Workspace/DynamicSys/pyenv_gpds/bin/activate"
export PYTHONPATH="$HOME/Workspace/DynamicSys/dynamic-systems/comparison/"

Rscript R/fn-model-comp.R
