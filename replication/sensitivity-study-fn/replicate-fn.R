## Run this script to reproduce FitzHugh-Nagumo (FN) system results in paper, for the same 100 random seeds

## Perform inference for the 100 seeds (can be run in parallel)
## Summarize results, make table and figures

source("fn-summarize.R")
library(parallel)
mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 4 ", i, " 21"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill4-nobs21/")

for (i in 1:100) {
  system(paste0("Rscript run-fn.R ", i), wait = TRUE)
}
library(parallel)
mclapply(1:2000, function(i){
  system("Rscript run-fn.R", wait=TRUE)
}, mc.cores=120)


mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 0 ", i, " 41"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill0-nobs41/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 1 ", i, " 41"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill1-nobs41/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 2 ", i, " 41"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill2-nobs41/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 3 ", i, " 41"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill3-nobs41/")

{
mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 4 ", i, " 11"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill4-nobs11/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 5 ", i, " 11"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill5-nobs11/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 6 ", i, " 11"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill6-nobs11/")

}

{
  mclapply(1:100, function(i){
    system(paste0("Rscript run-fn.R 4 ", i, " 21"), wait=TRUE)
  }, mc.cores=100)
  summarize("../results/fn-fill4-nobs21/")
  
  mclapply(1:100, function(i){
    system(paste0("Rscript run-fn.R 5 ", i, " 21"), wait=TRUE)
  }, mc.cores=100)
  summarize("../results/fn-fill5-nobs21/")
  
  mclapply(1:100, function(i){
    system(paste0("Rscript run-fn.R 3 ", i, " 21"), wait=TRUE)
  }, mc.cores=100)
  summarize("../results/fn-fill3-nobs21/")
}

x <- mclapply(1:100, function(i){
  system(paste0("python3 fnmap.py --nobs 6 --fill 5 --seed_id ", i), wait=TRUE)
}, mc.cores=100)

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 5 ", i, " 6"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill5-nobs6/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 6 ", i, " 6"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill6-nobs6/")

mclapply(1:100, function(i){
  system(paste0("Rscript run-fn.R 7 ", i, " 6"), wait=TRUE)
}, mc.cores=100)
summarize("../results/fn-fill7-nobs6/")

x <- mclapply(1:100, function(i){
  system(paste0("python3 mle_euler.py --nobs 6 --fill 5 --seed_id ", i), wait=TRUE)
}, mc.cores=100)

x <- mclapply(1:100, function(i){
  system(paste0("python3 mle_fdgm.py --nobs 6 --fill 4 --seed_id ", i), wait=TRUE)
}, mc.cores=25)
