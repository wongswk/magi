library(gpds)

nobs.candidates <- c(5, 11, 26, 51, 101, 201, 401)
noise.candidates <- c(0.01, 0.1, 0.2, 0.5, 1.0, 2) * 5
filllevel.candidates <- 0:4

indicatorArray <- array(FALSE, dim=c(length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1 # --array=0-9999

indicatorArray[arg] <- TRUE

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 2, any)],
  noise = noise.candidates[apply(indicatorArray, 1, any)],
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 200,
  n.iter = 10000,
  burninRatio = 0.2,
  stepSizeFactor = 0.1,
  filllevel = filllevel.candidates[apply(indicatorArray, 3, any)],
  modelName = "hes1"
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(config$ndis <= 801){
  for(dummy in 1:20){
    config$seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9
    source("R/hes1-model.R")    
  }
}
