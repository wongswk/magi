library(gpds)

nobs.candidates <- c(401)
filllevel.candidates <- 0
temperPrior.candidates <- c(TRUE)

indicatorArray <- array(FALSE, dim=c(length(temperPrior.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1

indicatorArray[arg] <- TRUE

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 2, any)],
  noise = c(0.15, 0.07),
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 1e4,
  burninRatio = 0.50,
  stepSizeFactor = 1,
  filllevel = filllevel.candidates[apply(indicatorArray, 3, any)],
  modelName = "FN",
  startXAtTruth = FALSE,
  startThetaAtTruth = FALSE,
  startSigmaAtTruth = FALSE,
  useGPmean = TRUE,
  forseTrueMean = FALSE,
  phase2 = FALSE,
  phase3 = FALSE,
  temperPrior = temperPrior.candidates[apply(indicatorArray, 1, any)],
  max.epoch = 9,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1]
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(config$ndis <= 201){
  source("R/fn-model.R")  
}
