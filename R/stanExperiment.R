library(gpds)

nobs.candidates <- c(51, 101, 201, 401, 801, 1601)
noise.candidates <- c(0.001, 0.5, 0.5, 0.5)
filllevel.candidates <- 0:5
finiteDifferenceType.candidates <- 0:1

indicatorArray <- array(FALSE, dim=c(length(finiteDifferenceType.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1
indicatorArray[arg] <- TRUE

seed <- (which(apply(indicatorArray, 1, any)) * 1 +
           which(apply(indicatorArray, 2, any)) * 10 +
           which(apply(indicatorArray, 3, any)) * 100)
  
config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = rep(noise.candidates[apply(indicatorArray, 2, any)], 2),
  seed = seed,
  npostplot = 50,
  filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
  modelName = "FN"
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

stanConfig <- list(
  sigma_obs=0.5,
  sigma_xdot=0.1,
  finiteDifferenceType=finiteDifferenceType.candidates[apply(indicatorArray, 1, any)]
)

if(config$ndis <= 1601){
  source("R/m-finiteDifference-stan-simu.R")
}
