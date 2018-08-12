#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 7 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 4000 # Memory request
#SBATCH -t 1-01:20 # (D-HH:MM)
#SBATCH -o /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.out # Standard output
#SBATCH -e /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.err # Standard error

Rscript R/stanExperiment.R ${SLURM_ARRAY_TASK_ID}
