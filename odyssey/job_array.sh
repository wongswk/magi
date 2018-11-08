#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 4000 # Memory request
#SBATCH -t 1-01:20 # (D-HH:MM)
#SBATCH -o /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.out # Standard output
#SBATCH -e /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.err # Standard error

# if (( $SLURM_ARRAY_TASK_ID % 4 == 0 ))
# then
#   Rscript R/hes1-log-model.R ${SLURM_ARRAY_TASK_ID}
# fi
# 
# if (( $SLURM_ARRAY_TASK_ID % 4 == 1 ))
# then
#   Rscript R/hes1-log-async-partial-observations.R ${SLURM_ARRAY_TASK_ID}
# fi
# 
# if (( $SLURM_ARRAY_TASK_ID % 4 == 2 ))
# then
#   Rscript R/hes1-model.R ${SLURM_ARRAY_TASK_ID}
# fi
# 
# if (( $SLURM_ARRAY_TASK_ID % 4 == 3 ))
# then
#   Rscript R/hes1-async-partial-observations.R ${SLURM_ARRAY_TASK_ID}
# fi

Rscript R/largeParameterExperiment.R ${SLURM_ARRAY_TASK_ID}
# Rscript R/fn-ramsay.R ${SLURM_ARRAY_TASK_ID}
