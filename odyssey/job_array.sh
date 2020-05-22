#!/bin/bash
#SBATCH -J gpdsFixFG # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 64 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 256000 # Memory request
#SBATCH -t 1-03:20 # (D-HH:MM)
#SBATCH -o /n/holyscratch01/kou_lab/shihaoyang/dynamic_sys/summary_%a.out # Standard output
#SBATCH -e /n/holyscratch01/kou_lab/shihaoyang/dynamic_sys/summary_%a.err # Standard error

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

# Rscript R/largeParameterExperiment.R ${SLURM_ARRAY_TASK_ID}
# Rscript R/fn-ramsay.R ${SLURM_ARRAY_TASK_ID}
# Rscript R/addRamsayFN.R ${SLURM_ARRAY_TASK_ID}
# Rscript R/hes1-log-async-partialobs-cpp-temper.R
Rscript R/hes1-log-7param-summary.R ${SLURM_ARRAY_TASK_ID}
