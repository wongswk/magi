## Run this script to reproduce FitzHugh-Nagumo (FN) system results in paper, for the same 100 random seeds

## Perform inference for the 100 seeds (can be run in parallel)
## Summarize results, make table and figures

source("fn-summarize.R")
for (i in 1:100) {
  system(paste0("Rscript run-fn.R 4 ", i, " 21"), wait=FALSE)
}
summarize("../results/fn-fill4-nobs21/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 0 ", i, " 41"), wait=FALSE)
}
summarize("../results/fn-fill0-nobs41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 1 ", i, " 41"), wait=FALSE)
}
summarize("../results/fn-fill1-nobs41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 2 ", i, " 41"), wait=FALSE)
}
summarize("../results/fn-fill2-nobs41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 3 ", i, " 41"), wait=FALSE)
}
summarize("../results/fn-fill3-nobs41/")

for (i in 1:100) {
  system(paste0("python fnmap.py --seed_id ", i), wait=FALSE)
}
