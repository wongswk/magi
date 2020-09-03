## Run this script to reproduce protein transduction results in paper, for the same 100 random seeds

## Perform inference for the 100 seeds (can be run in parallel)
for (i in 1:100) {
  system(paste0("R --vanilla --args ", i, " --no-save < run-ptrans-highnoise.R"))
}

## Summarize results, make table and figures
source("ptrans-summarize-highnoise.R")