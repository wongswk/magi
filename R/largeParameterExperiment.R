library(gpds)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2)
filllevel.candidates <- 0:6
temperPrior.candidates <- c(TRUE, FALSE)

indicatorArray <- array(FALSE, dim=c(length(temperPrior.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1
arg <- 67
indicatorArray[arg] <- TRUE

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = rep(noise.candidates[apply(indicatorArray, 2, any)], 2),
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 5000,
  burninRatio = 0.20,
  stepSizeFactor = 1,
  filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  forseTrueMean = TRUE,
  phase2 = FALSE,
  phase3 = FALSE,
  temperPrior = temperPrior.candidates[apply(indicatorArray, 1, any)],
  max.epoch = 0,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(config$ndis <= 401){
  source("R/fn-model.R")  
}
