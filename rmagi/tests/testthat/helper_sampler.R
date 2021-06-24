#' chain Sampler generate one MCMC chain
#' 
#' @param config need at least `n.iter`, `burninRatio`, `loglikflag`, `hmcSteps`
#' 
#' @export 
chainSampler <- function(config, xInit, singleSampler, stepLowInit, verbose=TRUE, 
                         thetaDim=3){
  numparam <- length(xInit)
  n.iter <- config$n.iter
  xth.formal <- matrix(NA, n.iter, numparam)
  stepLow <- stepLowInit
  accepts <- vector("numeric", n.iter)
  accepts[1] <- 0
  lliklist <-  vector("numeric", n.iter)
  xth.formal[1, ] <- xInit
  burnin <- as.integer(n.iter*config$burninRatio)
  
  for (t in 2:n.iter) {
    rstep <- runif(length(stepLow), stepLow, 2*stepLow)
    foo <- singleSampler(xth.formal[t-1,], rstep)
    
    xth.formal[t,] <- foo$final
    accepts[t] <- foo$acc
    
    if (t < burnin & t > 10) {
      if (mean(tail(accepts[1:t],100)) > 0.9) {
        stepLow <- stepLow * 1.005
      } else if (mean(tail(accepts[1:t],100)) < 0.6) {
        stepLow <- stepLow * .995
      }
      if( t %% 100 ==0) {
        xthsd <- apply(tail(xth.formal[1:t,],100), 2, sd)
        if(mean(xthsd)>0) stepLow <- 0.05*xthsd/mean(xthsd)*mean(stepLow) + 0.95*stepLow
      }
    }
    lliklist[t] <- foo$lpr
    
    if(verbose && t %% 100 ==0) 
      methods::show(c(t, mean(tail(accepts[1:t],100)), tail(foo$final, thetaDim)))
    
  }
  list(xth=xth.formal, lliklist=lliklist, stepLow=stepLow)
}

#' pilot Sampler generate one MCMC chain
#' 
#' @param config need at least `n.iter`, `burninRatio`, 
#' `loglikflag`, `hmcSteps`, `pilotRatio`, `stepSize`, 
#' 
#' @export 
runPilot <- function(pilotIndex, xsim.obs, r.nobs, signr.nobs, curphi,
                     gpsmoothFuncList, thetaSize, config){
  xsim.obs <- xsim.obs[pilotIndex, ]
  r.nobs <- r.nobs[pilotIndex, pilotIndex]
  signr.nobs <- signr.nobs[pilotIndex, pilotIndex]
  
  pilotCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r.nobs, signr.nobs, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  
  pilotIterations <- as.integer(config$n.iter*config$pilotRatio)
  config$n.iter <- pilotIterations
  burnin <- as.integer(config$n.iter*config$burninRatio)
  
  xInit <- lapply(gpsmoothFuncList, function(f) f(xsim.obs$time))
  xInit <- unlist(xInit)
  xInit <- c(xInit, rep(1, thetaSize))
  
  stepLowInit <- rep(config$stepSize, length(data.matrix(xsim.obs[,-1])) + thetaSize)
  
  singleSampler <- function(xthetaValues, stepSize) 
    xthetaSample(data.matrix(xsim.obs[,-1]), pilotCov, cursigma, 
                 xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
                 priorTemperature = config$priorTemperature, modelName = config$modelName)
  # save(list=ls(), file="debug_in_fun.rda")
  chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)
  
  gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim.obs[,-1]))+1):(ncol(chainSamplesOut$xth))],
                xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim.obs[,-1]))], 
                               dim=c(config$n.iter-burnin, nrow(xsim.obs), ncol(xsim.obs)-1)))
  list(theta = colMeans(gpode$theta), x = apply(gpode$xsampled, 2:3, mean))
}
