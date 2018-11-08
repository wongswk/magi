library(gpds)

config <- list(
  nobs = 401,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 1e4,
  burninRatio = 0.50,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = FALSE,
  startThetaAtTruth = FALSE,
  startSigmaAtTruth = FALSE,
  useGPmean = TRUE,
  forseTrueMean = FALSE,
  phase2 = FALSE,
  phase3 = FALSE,
  temperPrior = TRUE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1]
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

source("R/fn-model.R")  

