# TODO: functionalize verything here
# phi/sigma fitting method -----------------------------------------------------
phiAllMethodsExpr <- quote({
  cursigma <- rep(NA, ncol(xsim)-1)
  curphi <- matrix(NA, 2, ncol(xsim)-1)
  
  # MarginalLikelihood 
  for(j in 1:(ncol(xsim)-1)){
    fn <- function(par) -phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                      r.nobs, config$kernel)$value
    gr <- function(par) -as.vector(phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                                r.nobs, config$kernel)$grad)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarginalLikelihood <- cursigma
  curphiMarginalLikelihood <- curphi
  
  # LoocvLlik 
  for(j in 1:(ncol(xsim)-1)){
    fn <- function(par) -phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), 
                                           r.nobs, config$kernel)$value
    gr <- function(par) -as.vector(phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), 
                                                     r.nobs, config$kernel)$grad)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaLoocvLlik <- cursigma
  curphiLoocvLlik <- curphi
  
  # LoocvMse 
  for(j in 1:(ncol(xsim)-1)){
    fn <- function(par) -phisigloocvmseC( c(par, pram.true$sigma[j]) , data.matrix(xsim.obs[,1+j]), 
                                          r.nobs, config$kernel)$value
    gr <- function(par) -as.vector(phisigloocvmseC( c(par, pram.true$sigma[j]), data.matrix(xsim.obs[,1+j]), 
                                                    r.nobs, config$kernel)$grad)[1:2]
    marlikmap <- optim(rep(100, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- pram.true$sigma[j]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaLoocvMse <- cursigma
  curphiLoocvMse <- curphi
  
  # marllik+loocvllik+fftprior 
  for(j in 1:(ncol(xsim)-1)){
    priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
    
    fn <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                       max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
      -((marlik$value + loocvlik$value)/2 + penalty)
    }
    gr <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      grad <- -as.vector(marlik$grad + loocvlik$grad)/2
      grad[2] <- grad[2] + (par[2] - max(xsim.obs$time)*priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
      grad
    }
    testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikLoocvllkFftprior <- cursigma
  curphiMarllikLoocvllkFftprior <- curphi
  
  # marllik+fftprior 
  for(j in 1:(ncol(xsim)-1)){
    priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
    
    fn <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                       max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
      -(marlik$value + penalty)
    }
    gr <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      grad <- -as.vector(marlik$grad)
      grad[2] <- grad[2] + (par[2] - max(xsim.obs$time)*priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
      grad
    }
    testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikFftprior <- cursigma
  curphiMarllikFftprior <- curphi
  
  # marllik+loocvllik+fftGamaprior 
  for(j in 1:(ncol(xsim)-1)){
    priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
    
    desiredMode <- priorFactor["meanFactor"]
    betaRate <- uniroot(function(betaRate) pgamma(1, 1 + desiredMode*betaRate, betaRate)-0.95,
                        c(1e-3, 1e3))$root
    alphaRate <- 1 + desiredMode*betaRate
    # plot.function(function(x) dgamma(x, alphaRate, betaRate/5), 0, 5, n=1e3)
    fn <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      penalty <- dgamma(par[2], alphaRate, betaRate/max(xsim.obs$time), log=TRUE)
      -((marlik$value + loocvlik$value)/2 + penalty)
    }
    gr <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      grad <- -as.vector(marlik$grad + loocvlik$grad)/2
      grad[2] <- grad[2] - ((alphaRate-1)/par[2] - betaRate/max(xsim.obs$time))
      grad
    }
    testthat::expect_equal(gr(c(5,20,1))[2], (fn(c(5,20+1e-8,1)) - fn(c(5,20,1)))/1e-8, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikLoocvllkFftGammaprior <- cursigma
  curphiMarllikLoocvllkFftGammaprior <- curphi
  
  # marllik+fftGammaprior 
  for(j in 1:(ncol(xsim)-1)){
    priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
    
    desiredMode <- priorFactor["meanFactor"]
    betaRate <- uniroot(function(betaRate) pgamma(1, 1 + desiredMode*betaRate, betaRate)-0.95,
                        c(1e-3, 1e3))$root
    alphaRate <- 1 + desiredMode*betaRate
    
    fn <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      penalty <- dgamma(par[2], alphaRate, betaRate/max(xsim.obs$time), log=TRUE)
      -(marlik$value + penalty)
    }
    gr <- function(par) {
      marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
      grad <- -as.vector(marlik$grad)
      grad[2] <- grad[2] - ((alphaRate-1)/par[2] - betaRate/max(xsim.obs$time))
      grad
    }
    testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikFftGammaprior <- cursigma
  curphiMarllikFftGammaprior <- curphi
  rm(cursigma, curphi)
})

# profile likelihoood for phi2 ---------------------------------------------------
checkPhi2FitMarllik <- function(phi2, j=3){
  fn <- function(par) {
    phisig <- c(par[1], phi2, par[2])
    marlik <- phisigllikC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -marlik$value
  }
  gr <- function(par) {
    phisig <- c(par[1], phi2, par[2])
    marlik <- phisigllikC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -as.vector(marlik$grad)[-2]
  }
  marlikmap <- optim(rep(100, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  par <- marlikmap$par
  c(-marlikmap$value, par[1], phi2, par[2])
}

checkPhi2FitLoocvllik <- function(phi2, j=3){
  fn <- function(par) {
    phisig <- c(par[1], phi2, par[2])
    marlik <- phisigloocvllikC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -marlik$value
  }
  gr <- function(par) {
    phisig <- c(par[1], phi2, par[2])
    marlik <- phisigloocvllikC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -as.vector(marlik$grad)[-2]
  }
  marlikmap <- optim(rep(100, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  par <- marlikmap$par
  c(-marlikmap$value, par[1], phi2, par[2])
}

checkPhi2FitLoocvmse <- function(phi2, j=3){
  fn <- function(par) {
    phisig <- c(par[1], phi2, pram.true$sigma[j])
    marlik <- phisigloocvmseC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -marlik$value
  }
  gr <- function(par) {
    phisig <- c(par[1], phi2, pram.true$sigma[j])
    marlik <- phisigloocvmseC( phisig, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    -as.vector(marlik$grad)[1]
  }
  marlikmap <- optim(rep(100, 1), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  par <- marlikmap$par
  c(-marlikmap$value, par[1], phi2, pram.true$sigma[j])
}

getGPsmoothFunc <- function(thisphi, thissigma, j, showplot=TRUE){
  ynew <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                       t(thisphi), t(thissigma), kerneltype=config$kernel)
  gpsmoothFunc <- approxfun(xsim$time, ynew)
  if(showplot){
    plot.function(gpsmoothFunc, from = min(xsim$time), to = max(xsim$time),
                  lty = 2, col = j, add = TRUE)  
  }
  gpsmoothFunc
}

# MLE optim functions ----------------------------------------------------------------
fulloptim <- function(xInit, thetaInit, curphi, cursigma){
  fn <- function(par) -llikXthetaphisigma( par )$value
  gr <- function(par) -as.vector(llikXthetaphisigma( par )$grad)
  marlikmap <- optim(c(xInit, thetaInit, curphi, cursigma), 
                     fn, gr, method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  
  xInit[] <- marlikmap$par[xId]
  thetaInit[] <- marlikmap$par[thetaId]
  curphi[] <- marlikmap$par[phiId]
  cursigma[] <- marlikmap$par[sigmaId]
  list(xInit = xInit, 
       thetaInit = thetaInit, 
       curphi = curphi, 
       cursigma = cursigma)
}

phisigmaoptim <- function(xInit, thetaInit, curphi, cursigma){
  fullInit <- c(xInit, thetaInit, curphi, cursigma)
  fn <- function(par) {
    fullInit[c(phiId, sigmaId)] <- par
    -llikXthetaphisigma( fullInit )$value
  }
  gr <- function(par) {
    fullInit[c(phiId, sigmaId)] <- par
    -as.vector(llikXthetaphisigma( fullInit )$grad[c(phiId, sigmaId)])
  }
  marlikmap <- optim(c(curphi, cursigma), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  curphi[] <- marlikmap$par[1:length(curphi)]
  cursigma[] <- marlikmap$par[-(1:length(curphi))]
  list(curphi = curphi,
       cursigma = cursigma)
}

thetaoptim <- function(xInit, thetaInit, curphi, cursigma){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    -xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "Hes1" )$value
  }
  gr <- function(par) {
    -as.vector(xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "Hes1" )$grad[-(1:length(xInit))])
  }
  marlikmap <- optim(c(thetaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  thetaInit[] <- marlikmap$par
  list(thetaInit = thetaInit)
}

xthetaoptim <- function(xInit, thetaInit, curphi, cursigma, priorTemperature = rep(1,2)){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    -xthetallikRcpp( yobs, curCov, cursigma, par, "Hes1", priorTemperatureInput=priorTemperature )$value
  }
  gr <- function(par) {
    -as.vector(xthetallikRcpp( yobs, curCov, cursigma, par, "Hes1", priorTemperatureInput=priorTemperature )$grad)
  }
  marlikmap <- optim(c(xInit, thetaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  xInit[] <- marlikmap$par[1:length(xInit)]
  thetaInit[] <- marlikmap$par[-(1:length(xInit))]
  list(xInit = xInit,
       thetaInit = thetaInit)
}

xthetasigmaoptim <- function(xInit, thetaInit, curphi, sigmaInit, priorTemperature = rep(1,2)){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    xthis <- matrix(par[xId], ncol=ncol(yobs))
    thetathis <- par[thetaId]
    sigmathis <- par[(max(thetaId)+1):length(par)]
    -xthetasigmallikRcpp( xthis, thetathis, sigmathis, yobs, curCov, priorTemperature, FALSE, "Hes1" )$value
  }
  gr <- function(par) {
    xthis <- matrix(par[xId], ncol=ncol(yobs))
    thetathis <- par[thetaId]
    sigmathis <- par[(max(thetaId)+1):length(par)]
    -as.vector(xthetasigmallikRcpp( xthis, thetathis, sigmathis, yobs, curCov, priorTemperature, FALSE, "Hes1" )$grad)
  }
  marlikmap <- optim(c(xInit, thetaInit, sigmaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  xInit[] <- marlikmap$par[xId]
  thetaInit[] <- marlikmap$par[thetaId]
  sigmaInit[] <- marlikmap$par[(max(thetaId)+1):length(marlikmap$par)]
  list(xInit = xInit,
       thetaInit = thetaInit,
       cursigma = sigmaInit)
}

x3thetaphi3optim <- function(xInit, thetaInit, curphi, cursigma){
  fullInit <- c(xInit, thetaInit, curphi, cursigma)
  
  x3Id <- (length(xInit[, -3]) + 1):length(xInit)
  thetaId <- (max(x3Id)+1):(max(x3Id)+length(thetaInit))
  phi3Id <- (max(thetaId) + length(curphi[,-3]) + 1):(max(thetaId) + length(curphi))
  
  fn <- function(par) {
    fullInit[c(x3Id, thetaId, phi3Id)] <- par
    -llikXthetaphisigma( fullInit )$value
  }
  gr <- function(par) {
    fullInit[c(x3Id, thetaId, phi3Id)] <- par
    -as.vector(llikXthetaphisigma( fullInit )$grad[c(x3Id, thetaId, phi3Id)])
  }
  marlikmap <- optim(c(xInit[, 3], thetaInit, curphi[, 3]), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  xInit[, 3] <- marlikmap$par[1:length(x3Id)]
  thetaInit <- marlikmap$par[(length(x3Id)+1):(length(x3Id)+length(thetaId))]
  curphi[, 3] <- tail(marlikmap$par, length(phi3Id))
  list(xInit = xInit,
       thetaInit = thetaInit,
       curphi = curphi)
}

llikXthetaphisigma <- function(xthetaphisigma) {
  xInitial <- matrix(xthetaphisigma[xId], nrow=obsDim[1], ncol=obsDim[2])
  thetaInitial <- xthetaphisigma[thetaId]
  phiInitial <- matrix(xthetaphisigma[phiId], nrow=2)
  sigmaInitial <- xthetaphisigma[sigmaId]
  xthetaphisigmallikRcpp(xInitial, thetaInitial, phiInitial, sigmaInitial,
                         yobs, xsim$time, modelName = "Hes1")
}

x3thetaphi3sgd <- function(xInit, thetaInit, curphi, cursigma, maxit=1e4, learningRateInput=1e-8){
  fullInit <- c(xInit, thetaInit, curphi, cursigma)
  
  x3Id <- (length(xInit[, -3]) + 1):length(xInit)
  thetaId <- (max(x3Id)+1):(max(x3Id)+length(thetaInit))
  phi3Id <- (max(thetaId) + length(curphi[,-3]) + 1):(max(thetaId) + length(curphi))
  
  fg <- function(par) {
    fullInit[c(x3Id, thetaId, phi3Id)] <- par
    llik <- llikXthetaphisigma( fullInit )
    list(fn=-llik$value, gr=-as.vector(llik$grad[c(x3Id, thetaId, phi3Id)]))
  }
  par <- c(xInit[, 3], thetaInit, curphi[, 3])
  fn.last <- fg(par)$fn
  par.last <- par
  for (i in 1:maxit){
    fgv <- fg(par)
    if(fgv$fn - fn.last > abs(fn.last) * 0.10){
      par <- par.last
      learningRateInput = learningRateInput / 10
      next
    }
    if(i %% 100 == 0){
      learningRateInput = learningRateInput * 10
    }
    grValue <- fgv$gr
    learningRate <- 0.001 * max(abs(tail(par, length(phi3Id)) / tail(grValue, length(phi3Id))))
    learningRate <- min(learningRate, learningRateInput)
    if(learningRate < 1e-20){
      break
    }
    par.last <- par
    fn.last <- fgv$fn
    par <- par - learningRate * grValue
    par <- pmax(par, 0.001)
    cat("------ ", fgv$fn)
    print(learningRate)
    print(round(par[(length(x3Id)+1):(length(x3Id)+length(thetaId))], 3))
    print(tail(par, length(phi3Id)))
  }
  xInit[, 3] <- par[1:length(x3Id)]
  thetaInit <- par[(length(x3Id)+1):(length(x3Id)+length(thetaId))]
  curphi[, 3] <- tail(par, length(phi3Id))
  list(xInit = xInit,
       thetaInit = thetaInit,
       curphi = curphi)
}
