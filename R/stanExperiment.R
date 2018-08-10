library(gpds)

nobs.candidates <- c(101, 201)
noise.candidates <- c(0.001, 0.5, 0.5, 0.5)
filllevel.candidates <- 0:3
maternDf.candidates <- seq(2.1, 4.5, 0.1)

indicatorArray <- array(FALSE, dim=c(length(maternDf.candidates), 
                                     length(noise.candidates), 
                                     length(nobs.candidates),
                                     length(filllevel.candidates) ))
arg <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(arg) %% length(indicatorArray) + 1
indicatorArray[arg] <- TRUE

seed <- (which(apply(indicatorArray, 2, any)) * 10 +
           which(apply(indicatorArray, 3, any)) * 100)
  
thisDf <- maternDf.candidates[apply(indicatorArray, 1, any)]

config <- list(
  nobs = nobs.candidates[apply(indicatorArray, 3, any)],
  noise = rep(noise.candidates[apply(indicatorArray, 2, any)], 2),
  seed = seed,
  npostplot = 50,
  filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
  modelName = "FN",
  kernel = paste0("generalMatern-", round(thisDf, 2))
)
config$ndis <- (config$nobs-1)*2^config$filllevel+1

stanConfig <- list(
  sigma_obs=0.5,
  sigma_xdot=0.1
)

if(config$ndis <= 801){
  source("R/m-finiteDifference-stan-simu.R")
}
