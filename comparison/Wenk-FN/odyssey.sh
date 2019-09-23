#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 26 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 8000 # Memory request
#SBATCH -t 1-01:20 # (D-HH:MM)
#SBATCH -o /n/scratchlfs/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.out # Standard output
#SBATCH -e /n/scratchlfs/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.err # Standard error


source "$HOME/Workspace/DynamicSys/pyenv_gpds/bin/activate"
export PYTHONPATH="$HOME/Workspace/DynamicSys/dynamic-systems/comparison/"
mkdir -p "$PYTHONPATH/Wenk-FN/Wenk19-FN-${SLURM_ARRAY_TASK_ID}" && cd "$PYTHONPATH/Wenk-FN/Wenk19-FN-${SLURM_ARRAY_TASK_ID}"

python3 $PYTHONPATH/FGPGM/mainFiles/FitzHughNagumo/createExperiments.py --slurm_array_task_id $SLURM_ARRAY_TASK_ID
python3 $PYTHONPATH/FGPGM/mainFiles/FitzHughNagumo/getHyperparams.py
python3 $PYTHONPATH/FGPGM/mainFiles/FitzHughNagumo/doFGPGM.py
