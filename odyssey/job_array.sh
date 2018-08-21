#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 5000 # Memory request
#SBATCH -t 2-01:20 # (D-HH:MM)
#SBATCH -o /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystemDropout%a.out # Standard output
#SBATCH -e /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystemDropout%a.err # Standard error

Rscript R/dropoutExperiment.R ${SLURM_ARRAY_TASK_ID}
