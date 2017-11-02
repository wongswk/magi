#!/bin/bash
#SBATCH -J DynamicSystem # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p kou # Partition
#SBATCH --mem 1000 # Memory request
#SBATCH -t 0-06:00 # (D-HH:MM)
#SBATCH -o /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.out # Standard output
#SBATCH -e /n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem%a.err # Standard error

Rscript R/two-phase-repeated-sample-pluginTruth.R ${SLURM_ARRAY_TASK_ID}
