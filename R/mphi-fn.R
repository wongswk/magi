# sample from $\| f - m X\|_I$ --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = TRUE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1e5),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# remove force $K$ to be identity matrix from 1 --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1e5),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# plug-in true curve as mean --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# plug-in biased phase8 inference curve as mean --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = FALSE,
  forsePhase8Mean = TRUE,
  priorTemperature = c(1, 1),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")


