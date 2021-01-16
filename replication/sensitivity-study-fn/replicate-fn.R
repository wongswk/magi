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

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-fixsigma-truemu.R 5 ", i, " 6 linear"), wait=FALSE)
}
summarize("../results/fn-sparse-fixsigma-truemu-fill5-nobs6-initXlinear/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparse-fixsigma-truemu.R 5 ", i, " 6 truth"), wait=FALSE)
}
summarize("../results/fn-sparse-fixsigma-truemu-fill5-nobs6-initXtruth/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-fixphi-fixsigma-truemu.R 5 ", i, " 6 linear"), wait=FALSE)
}
summarize("../results/fn-fixphi-fixsigma-truemu-fill5-nobs6-initXlinear/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-fixphi-fixsigma-truemu.R 5 ", i, " 6 truth"), wait=FALSE)
}
summarize("../results/fn-fixphi-fixsigma-truemu-fill5-nobs6-initXtruth/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-fill5-nobs6/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-init-truth.R 5 ", i, " 6"), wait=FALSE)
}
summarize("../results/fn-init-truth-fill5-nobs6/")

for(filllevel in 0:4){
  for (i in 1:4) {
    system(paste0("Rscript run-linear-fixphi-fixsigma.R ", filllevel, " ", i, " 41"), wait=FALSE)
  }
}

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 40"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-timeend40/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 41 40 161"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs41-ninterp161-timeend40/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 21 20 81"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21-ninterp81-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 21 20 161"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21-ninterp161-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 21 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 3 ", i, " 21 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill3-nobs21-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-noband-fixphi.R 3 ", i, " 21 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-fill3-nobs21-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 21 20 21"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21-ninterp21-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 21 20 33"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs21-ninterp33-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 21 20 21"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs21-ninterp21-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 20 161"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-ninterp161-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 41 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs41-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 21 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs21-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 21 20 81"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs21-ninterp81-timeend20/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 2 ", i, " 41 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill2-nobs41-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 3 ", i, " 41 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill3-nobs41-withmeanBand-fixphi-heating/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 4 ", i, " 11 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill4-nobs11-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 5 ", i, " 11 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill5-nobs11-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 6 ", i, " 11 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill6-nobs11-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 4 ", i, " 21 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill4-nobs21-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-ideal.R 3 ", i, " 21 heating withmeanBand"), wait=FALSE)
}
summarize("../results/fn-ideal-fill3-nobs21-withmeanBand-fixphi-heating/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 5 ", i, " 11 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill5-nobs11-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 11 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs11-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 6 ", i, " 11 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill6-nobs11-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 5 ", i, " 11 20 81"), wait=FALSE)
}
summarize("../results/fn-sparse-fill5-nobs11-ninterp81-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 11 20 81"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs11-ninterp81-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 6 ", i, " 11 20 81"), wait=FALSE)
}
summarize("../results/fn-sparse-fill6-nobs11-ninterp81-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 7 ", i, " 11 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill7-nobs11-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-interp.R 4 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill4-nobs11-ninterp41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-interp.R 5 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill5-nobs11-ninterp41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-interp.R 6 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill6-nobs11-ninterp41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-interp.R 4 ", i, " 11 81"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill4-nobs11-ninterp81/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 4 ", i, " 11 81"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill4-nobs11-reEstPhi-ninterp81/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 5 ", i, " 11 81"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill5-nobs11-reEstPhi-ninterp81/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 6 ", i, " 11 81"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill6-nobs11-reEstPhi-ninterp81/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 4 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill4-nobs11-reEstPhi-ninterp41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 5 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill5-nobs11-reEstPhi-ninterp41/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-sparseComponent-opt-reestphi-interp.R 6 ", i, " 11 41"), wait=FALSE)
}
summarize("../results/fn-sparseComponent-opt-interp-fill6-nobs11-reEstPhi-ninterp41/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 4 ", i, " 21 20 321"), wait=FALSE)
}
summarize("../results/fn-sparse-fill4-nobs21-ninterp321-timeend20/")


for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 20 161"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-ninterp161-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 2 ", i, " 41 20 41"), wait=FALSE)
}
summarize("../results/fn-sparse-fill2-nobs41-ninterp41-timeend20/")

for (i in 1:100) {
  system(paste0("Rscript run-fn-horizon.R 3 ", i, " 41 20 321"), wait=FALSE)
}
summarize("../results/fn-sparse-fill3-nobs41-ninterp321-timeend20/")

## Summarize results, make table and figures
# source("fn-summarize.R")
# summarize("~/Workspace/DynamicSys/results/sensitivity-result/fn-fill3-nobs21/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R 2 ", i, " 41"), wait=FALSE)
}
summarize("../results/fn-fill2-nobs41/")


for(seed_i in 1:5){
  for (i in 0:10) {
    system(paste0("Rscript run-fn.R ", i, " ", seed_i, " 41"), wait=FALSE)
  }
}
