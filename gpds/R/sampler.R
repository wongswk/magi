#' chain Sampler generate one MCMC chain
#' 
#' @param config need at least `n.iter`, `burninRatio`, `loglikflag`, `hmcSteps`
#' 
#' @export 
chainSampler <- function(config, xInit, singleSampler, stepLowInit, verbose=TRUE){
  numparam <- length(xInit)
  n.iter <- config$n.iter
  stepLow <- stepLowInit
  xth.formal <- matrix(NA, n.iter, numparam)
  accepts <- c(1)
  lliklist <- c()
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
      xthsd <- apply(tail(xth.formal[1:t,],n.iter*config$burninRatio*0.1), 2, sd)
      if(mean(xthsd)>0) stepLow <- 0.01*xthsd/mean(xthsd)*mean(stepLow) + 0.99*stepLow
    }
    lliklist[t] <- foo$lpr
    
    if(verbose && t %% 100 ==0) show(c(t, mean(tail(accepts[1:t],100)), tail(foo$final, 3)))
  }
  list(xth.formal, lliklist)
}