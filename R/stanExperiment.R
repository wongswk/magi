library(gpds)

nobs.candidates <- c(101, 201)
noise.candidates <- c(0.001, 0.5, 0.5, 0.5)
filllevel.candidates <- 0:3
kernel.candidates <- c("generalMatern-2.01", "generalMatern-2.05", "rbf", 
                       "finiteDifference1h", "finiteDifference2h")
kernel.maternDf <- paste0("generalMatern-", round(seq(2.1, 4.5, 0.1), 2))
kernel.candidates <- c(kernel.candidates, kernel.maternDf)

for(arg in 1:960){
  indicatorArray <- array(FALSE, dim=c(length(kernel.candidates), 
                                       length(noise.candidates), 
                                       length(nobs.candidates),
                                       length(filllevel.candidates) ))
  
  indicatorArray[arg] <- TRUE
  
  seed <- (which(apply(indicatorArray, 2, any)) * 10 +
             which(apply(indicatorArray, 3, any)) * 100)
  
  config <- list(
    nobs = nobs.candidates[apply(indicatorArray, 3, any)],
    noise = rep(noise.candidates[apply(indicatorArray, 2, any)], 2),
    seed = seed,
    npostplot = 50,
    filllevel = filllevel.candidates[apply(indicatorArray, 4, any)],
    modelName = "FN",
    kernel = kernel.candidates[apply(indicatorArray, 1, any)]
  )
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  
  stanConfig <- list(
    sigma_obs=0.5
  )
  
  if(config$ndis <= 801){
    source("R/m-finiteDifference-stan-simu.R")
  }
  
}
