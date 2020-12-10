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

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 5 ", i, " 6 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill5-nobs6-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 5 ", i, " 6 cooling withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill5-nobs6-withmeanBand-fixphi-cooling/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 5 ", i, " 6 notempering withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill5-nobs6-withmeanBand-fixphi-notempering/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse.R 3 ", i, " 21"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse.R 4 ", i, " 11"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs11/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-sparse-fill5-nobs6/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-cooling.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-sparse-cooling-fill5-nobs6/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-fixsigma.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-sparse-fixsigma-fill5-nobs6/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-fixsigma-cooling.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-sparse-fixsigma-cooling-fill5-nobs6/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-fixsigma-truemu.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-sparse-fixsigma-truemu-fill5-nobs6/")

## Summarize results, make table and figures
source("fn-summarize.R")
summarize("~/Workspace/DynamicSys/results/sensitivity-result/fn-fill3-nobs21/")

