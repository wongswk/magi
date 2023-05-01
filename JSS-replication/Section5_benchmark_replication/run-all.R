# Run the three methods on the Hes1 model, 100 repetitions each
# For simplicity this script runs everything sequentially on a single CPU core, but the repetitions can be run in parallel to reduce the time needed to complete the benchmark

nseeds <- 100

for (i in 1:nseeds) {
  system(paste0("Rscript collocInfer-hes1.R ", i))
}

for (i in 1:nseeds) {
  system(paste0("Rscript MAGI-hes1.R ", i))
}

for (i in 1:nseeds) {
  system(paste0("Rscript deBInfer-hes1.R ", i))
}


# Summarize results of the three methods
rmse_CollocInfer <- matrix(NA, nrow = nseeds, ncol = 3)
theta_CollocInfer <- matrix(NA, nrow = nseeds, ncol = 7)
timeUsedvec <- c()
for (seed in 1:nseeds) {
  if (file.exists(paste0("results/CollocInfer/Hes1-log-",seed,".rda"))) {
    load(file = paste0("results/CollocInfer/Hes1-log-",seed,".rda")) 
    
    rmse_CollocInfer[seed,] <- rmse_orig
    theta_CollocInfer[seed,] <- best.pars
    timeUsedvec[seed] <- TimeUsed
  }
}
apply(rmse_CollocInfer, 2, mean)
mean(timeUsedvec)/60

rmse_MAGI <- matrix(NA, nrow = nseeds, ncol = 3)
theta_MAGI <- matrix(NA, nrow = nseeds, ncol = 7)
timeUsedvec <- c()

for (seed in 1:nseeds) {
  if (file.exists(paste0("results/magi/Hes1-log-",seed,".rda"))) {
    load(file = paste0("results/magi/Hes1-log-",seed,".rda")) 
    
    rmse_MAGI[seed,] <- rmse_orig
    theta_MAGI[seed,] <- theta.est
    timeUsedvec[seed] <- TimeUsed
  }
}
apply(rmse_MAGI, 2, mean)
mean(timeUsedvec)/60

rmse_deBInfer <- matrix(NA, nrow = nseeds, ncol = 3)
theta_deBInfer <- matrix(NA, nrow = nseeds, ncol = 7)
timeUsedvec <- c()

for (seed in 1:nseeds) {
  if (file.exists(paste0("results/deBInfer/Hes1-log-",seed,".rda"))) {
    load(file = paste0("results/deBInfer/Hes1-log-",seed,".rda")) 
    
    rmse_deBInfer[seed,] <- rmse_orig
    theta_deBInfer[seed,] <- theta.est
    timeUsedvec[seed] <- TimeUsed
  }
}

apply(rmse_deBInfer, 2, mean)
mean(timeUsedvec)/60


