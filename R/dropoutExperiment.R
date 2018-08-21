library(gpds)

beta.candidates <- c(seq(1,13,2), Inf)
kernel.candidates <- c("finiteDifference1h", "finiteDifference2h", "rbf")
kernel.maternDf <- paste0("generalMatern-", round(c(2.01, 2.5, 3.5), 2))
kernel.candidates <- c(kernel.candidates, kernel.maternDf)
xCxTemperature.candidates <- c(1e6, 1)

indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                     length(xCxTemperature.candidates) ))

arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1
indicatorArray[arg] <- TRUE

seed <- 100

config <- list(
  nobs = 101,
  noise = c(0.5, 0.5),
  overrideNoise = TRUE,
  kernel = kernel.candidates[apply(indicatorArray, 1, any)],
  forceDiagKphi = FALSE,
  forceMean = c("gpsmooth", "truth", "phase8", "zero")[4],
  priorTemperature = c(1, xCxTemperature.candidates[apply(indicatorArray, 2, any)]),
  seed = 100,
  npostplot = 20,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 5000,
  burninRatio = 0.20,
  stepSizeFactor = 1,
  filllevel = 2,
  dropoutRate = 0.5,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  sigma_xdot = 0.1
)

config$ndis <- (config$nobs-1)*2^config$filllevel+1

source("R/finiteDifference2h-dropout.R")
