library(gpds)

nobs.candidates <- c(101, 201)
forceDiagKphi.candidates <- c(FALSE, TRUE)
filllevel.candidates <- 0:3
xCxTemperature.candidates <- c(1e6, 1)
kernel.candidates <- c("finiteDifference1h", "finiteDifference2h", "rbf")
forceMean.candidates <- c("gpsmooth", "truth", "phase8", "zero")
kernel.maternDf <- paste0("generalMatern-", round(c(2.01, 2.05, 2.1, seq(2.25, 4.5, 0.25)), 2))
kernel.candidates <- c(kernel.candidates, kernel.maternDf)

indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                     length(forceDiagKphi.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates),
                                     length(forceMean.candidates),
                                     length(xCxTemperature.candidates) ))

arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1
indicatorArray[arg] <- TRUE

seed <- which(apply(indicatorArray, 3, any)) * 100

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = c(0.5, 0.5),
  overrideNoise = TRUE,
  kernel = kernel.candidates[apply(indicatorArray, 1, any)],
  forceDiagKphi = forceDiagKphi.candidates[apply(indicatorArray, 2, any)],
  forceMean = forceMean.candidates[apply(indicatorArray, 5, any)],
  priorTemperature = c(1, xCxTemperature.candidates[apply(indicatorArray, 6, any)]),
  seed = seed,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 1000,
  n.iter = 5000,
  burninRatio = 0.20,
  stepSizeFactor = 1,
  filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  sigma_xdot = 0.1
)

config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  if(cofig$forceDiagKphi == FALSE){
    stop("finiteDifference but forceDiagKphi is FALSE")
  }
  if(forceMean == "gpsmooth"){
    stop("finiteDifference but forceMean is gpsmooth")
  }
}
  
if(config$ndis <= 801){
  source("R/fn-bias-check.R")
}
