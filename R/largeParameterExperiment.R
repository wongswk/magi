library(gpds)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:4
kernel.candidates <- c("generalMatern", "matern") # basically df in matern
nrep <- 100

indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %/% nrep + 1

indicatorArray[arg] <- TRUE

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = noise.candidates[apply(indicatorArray, 2, any)],
  kernel = kernel.candidates[apply(indicatorArray, 1, any)],
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 5,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 200,
  n.iter = 10000,
  burninRatio = 0.1,
  stepSizeFactor = 0.1,
  filllevel = filllevel.candidates[apply(indicatorArray, 4, any)]
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(config$ndis <= 801){
  source("R/two-phase-repeated-sample-pluginPostSample.R")  
}
