## Run this script to reproduce FitzHugh-Nagumo (FN) system results in paper, for the same 100 random seeds

## Perform inference for the 100 seeds (can be run in parallel)
for (i in 1:112) {
  system(paste0("Rscript run-fn.R ", i))
}

## Summarize results, make table and figures
source("fn-summarize.R")
