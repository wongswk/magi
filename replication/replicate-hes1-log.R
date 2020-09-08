## Run this script to reproduce Hes1 oscillation system results in paper, for the same 2000 random seeds

## Perform inference for the 2000 seeds (can be run in parallel)
for (i in 1:2000) {
  system(paste0("Rscript run-hes1-log.R ", i))
}

## Summarize results, make table and figures
source("hes1-log-summarize.R")
