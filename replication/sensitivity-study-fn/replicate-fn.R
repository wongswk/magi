## Run this script to reproduce FitzHugh-Nagumo (FN) system results in paper, for the same 100 random seeds

## Perform inference for the 100 seeds (can be run in parallel)
for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 5 ", i, " 6 cooling withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill5-nobs6-withmeanBand-fixphi-cooling/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 5 ", i, " 6 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill5-nobs6-withmeanBand-fixphi-heating/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 6 ", i, " 6 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill6-nobs6-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 6 ", i, " 6 cooling withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill6-nobs6-withmeanBand-fixphi-cooling/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 7 ", i, " 6 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill7-nobs6-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 7 ", i, " 6 cooling withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill7-nobs6-withmeanBand-fixphi-cooling/")

## Summarize results, make table and figures
source("fn-summarize.R")

