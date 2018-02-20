# TODO: functionalize verything here
phiAll3methodsExpr <- quote({
  cursigma <- rep(NA, ncol(xsim)-1)
  curphi <- matrix(NA, 2, ncol(xsim)-1)
  
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
})


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
