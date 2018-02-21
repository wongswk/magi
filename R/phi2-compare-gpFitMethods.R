#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
phi2compareMethods <- list()
configAll <- list()
nobs <- 11
mcresult <- parallel::mclapply(1:1000, function(dummyIterator){
  config <- list(
    nobs = nobs,
    noise = c(4,1,8)*0.2,
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 1e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = 3,
    modelName = "hes1"
  )
  
  
  config$ndis <- (config$nobs-1)*2^config$filllevel+1
  config$priorTemperature <- config$ndis / config$nobs
  if(grepl("/n/",getwd())){
    baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
  }else{
    baseDir <- "~/Workspace/DynamicSys/results/"  
  }
  outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                                "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                                "-ndis",ndis,"/"))
  system(paste("mkdir -p", outDir))
  
  pram.true <- list(
    theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
    x0 = c(0.5, 2, 1),
    phi = c(122.4027613, 41.8511396,  
            56.5612956, 91.4189948,
            164.3556832, 11.9474091),
    sigma = config$noise
  )
  times <- seq(0, 60*4, by = 0.01)
  
  modelODE <- function(t, state, parameters) {
    list(as.vector(gpds:::hes1modelODE(parameters, t(state))))
  }
  
  xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
  xtrue <- data.frame(xtrue)
  # if(dummyIterator < 4) matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
  
  xtrueFunc <- lapply(2:ncol(xtrue), function(j)
    approxfun(xtrue[, "time"], xtrue[, j]))
  
  xsim <- xtrue
  
  set.seed(config$seed)
  for(j in 1:(ncol(xsim)-1)){
    xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
  }
  
  xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
  # if(dummyIterator < 4) matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
  
  # if(dummyIterator < 4) matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)
  
  xsim <- insertNaN(xsim.obs,config$filllevel)
  
  tvec.full <- xsim$time
  tvec.nobs <- xsim.obs$time
  
  foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  signr <- -sign(foo)
  
  foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
  r.nobs <- abs(foo)
  r2.nobs <- r.nobs^2
  signr.nobs <- -sign(foo)
  
  cursigma <- rep(NA, ncol(xsim)-1)
  curphi <- matrix(NA, 2, ncol(xsim)-1)
  
  # MarginalLikelihood ---------------------------------------------------------
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
  cursigmaLoocvMarginalLikelihood <- cursigma
  curphiMarginalLikelihood <- curphi
  
  # LoocvLlik -------------------------------------------------------------------
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
  cursigmaLoocvLoocvLlik <- cursigma
  curphiLoocvLlik <- curphi
  
  # LoocvMse -------------------------------------------------------------------
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
  cursigmaLoocvLoocvMse <- cursigma
  curphiLoocvMse <- curphi
  
  # marllik+loocvllik+fftprior -------------------------------------------------
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
    # testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikLoocvllkFftprior <- cursigma
  curphiMarllikLoocvllkFftprior <- curphi
  
  # marllik+fftprior -----------------------------------------------------------
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
    # testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikFftprior <- cursigma
  curphiMarllikFftprior <- curphi
  
  # marllik+loocvllik+fftGamaprior -------------------------------------------------
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
    # testthat::expect_equal(gr(c(5,20,1))[2], (fn(c(5,20+1e-8,1)) - fn(c(5,20,1)))/1e-8, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikLoocvllkFftGammaprior <- cursigma
  curphiMarllikLoocvllkFftGammaprior <- curphi
  
  # marllik+fftGammaprior -----------------------------------------------------------
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
    # testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-4)
    marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaMarllikFftGammaprior <- cursigma
  curphiMarllikFftGammaprior <- curphi
  
  cursigmaLoocvLoocvMse
  cursigmaLoocvLoocvLlik
  cursigmaLoocvMarginalLikelihood
  
  curphiLoocvMse
  curphiLoocvLlik
  curphiMarginalLikelihood
  
  phi2compareMethods[[dummyIterator]] <-
    rbind(curphiLoocvMse[2,], curphiLoocvLlik[2,], curphiMarginalLikelihood[2,],
          curphiMarllikLoocvllkFftprior[2,], curphiMarllikFftprior[2,],
          curphiMarllikLoocvllkFftGammaprior[2,], curphiMarllikFftGammaprior[2,])
  
  configAll[[dummyIterator]] <- config    
  return(list(
    phi2compareMethods[[dummyIterator]],
    configAll[[dummyIterator]]
  ))
}, mc.cores = 8)

phi2compareMethods <- lapply(mcresult, function(x) x[[1]])
configAll <- lapply(mcresult, function(x) x[[2]])

phi2compareMethodsCube <- sapply(phi2compareMethods, identity, simplify = "array")
dimnames(phi2compareMethodsCube)[[1]] <- c("phi2LoocvMse", "phi2LoocvLlik", "phi2MarginalLikelihood",
                                           "marllik+loocvllik+fftNormalprior", "marllik+fftNormalprior",
                                           "marllik+loocvllik+fftGammaprior", "marllik+fftGammaprior")
dimnames(phi2compareMethodsCube)[[2]] <- c("hes1P", "hes1M", "hes1H")

ggobjs <- list()
for(component in c("hes1P", "hes1M", "hes1H")){
  x <- data.frame(log(t(phi2compareMethodsCube[,component,])))
  library(ggplot2);library(reshape2)
  data <- melt(x)
  ggobjs[[component]] <- 
    ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + 
    ggplot2::labs(title=component) + xlab("log phi2")
}

if(interactive()) print(ggobjs)

ml <- gridExtra::marrangeGrob(ggobjs, nrow=3, ncol=1, top="phi2")
ggsave(paste0("phi2-compare-gpFitMethods-nobs",nobs,".pdf"), ml, width=8, height = 16)
save.image(paste0("phi2-compare-gpFitMethods-nobs",nobs,".rda"))


# analysis after saving the data -----------------------------------------------
if(interactive()){
  nobs <- 26
  load(paste0("../results/2018-02-21/phi2-compare-gpFitMethods-nobs",nobs,".rda"))
  ggobjs <- list()
  phi2method <- c("phi2MarginalLikelihood", "marllik+loocvllik+fftNormalprior", "marllik+fftNormalprior")
  for(component in c("hes1P", "hes1M", "hes1H")){
    x <- data.frame(log(t(phi2compareMethodsCube[phi2method,component,])))
    library(ggplot2);library(reshape2)
    data <- melt(x)
    ggobjs[[component]] <- 
      ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + 
      ggplot2::labs(title=component) + xlab("log phi2")
  }
  ggobjs
  ml <- gridExtra::marrangeGrob(ggobjs, nrow=3, ncol=1, top="phi2")
  ggsave(paste0("../results/2018-02-21/phi2-compare-gpFitMethods-normalprior-nobs",nobs,".pdf"), ml, width=8, height = 16)
  
  ggobjs <- list()
  phi2method <- c("phi2MarginalLikelihood", "marllik+loocvllik+fftGammaprior", "marllik+fftGammaprior")
  for(component in c("hes1P", "hes1M", "hes1H")){
    x <- data.frame(log(t(phi2compareMethodsCube[phi2method,component,])))
    library(ggplot2);library(reshape2)
    data <- melt(x)
    ggobjs[[component]] <- 
      ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + 
      ggplot2::labs(title=component) + xlab("log phi2")
  }
  ggobjs
  ml <- gridExtra::marrangeGrob(ggobjs, nrow=3, ncol=1, top="phi2")
  ggsave(paste0("../results/2018-02-21/phi2-compare-gpFitMethods-gammaprior-nobs",nobs,".pdf"), ml, width=8, height = 16)
}
