#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
phi2compareMethods <- list()
configAll <- list()
for(dummyIterator in 1:1e3){
  config <- list(
    nobs = 11,
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
  matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
  
  xtrueFunc <- lapply(2:ncol(xtrue), function(j)
    approxfun(xtrue[, "time"], xtrue[, j]))
  
  xsim <- xtrue
  
  set.seed(config$seed)
  for(j in 1:(ncol(xsim)-1)){
    xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
  }
  
  xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
  
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)
  
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
  
  for(j in 1:(ncol(xsim)-1)){
    fn <- function(par) -phisigloocvllikC( c(par, pram.true$sigma[j]) , data.matrix(xsim.obs[,1+j]), 
                                           r.nobs, config$kernel)$value
    gr <- function(par) -as.vector(phisigloocvllikC( c(par, pram.true$sigma[j]), data.matrix(xsim.obs[,1+j]), 
                                                     r.nobs, config$kernel)$grad)[1:2]
    marlikmap <- optim(rep(100, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, 60*4*2, Inf))
    
    cursigma[j] <- pram.true$sigma[j]
    curphi[,j] <- marlikmap$par[1:2]
  }
  cursigmaLoocvLoocvMse <- cursigma
  curphiLoocvMse <- curphi
  
  cursigmaLoocvLoocvMse
  cursigmaLoocvLoocvLlik
  cursigmaLoocvMarginalLikelihood
  
  curphiLoocvMse
  curphiLoocvLlik
  curphiMarginalLikelihood
  
  phi2compareMethods[[dummyIterator]] <-
    rbind(curphiLoocvMse[2,], curphiLoocvLlik[2,], curphiMarginalLikelihood[2,])
  
  configAll[[dummyIterator]] <- config    
}

phi2compareMethodsCube <- sapply(phi2compareMethods, identity, simplify = "array")
dimnames(phi2compareMethodsCube)[[1]] <- c("phi2LoocvMse", "phi2LoocvLlik", "phi2MarginalLikelihood")
dimnames(phi2compareMethodsCube)[[2]] <- c("hes1P", "hes1M", "hes1H")


hist(phi2compareMethodsCube["phi2LoocvMse","hes1P",], breaks = 20, 
     col=rgb(1,0,0,0.3))
hist(phi2compareMethodsCube["phi2LoocvLlik","hes1P",], breaks = 20, 
     col=rgb(0,1,0,0.3), add=TRUE)
hist(phi2compareMethodsCube["phi2MarginalLikelihood","hes1P",], breaks = 20, 
     col=rgb(0,0,1,0.3), add=TRUE)

hist(phi2compareMethodsCube["phi2LoocvMse","hes1M",], breaks = 20, 
     col=rgb(1,0,0,0.3))
hist(phi2compareMethodsCube["phi2LoocvLlik","hes1M",], breaks = 20, 
     col=rgb(0,1,0,0.3), add=TRUE)
hist(phi2compareMethodsCube["phi2MarginalLikelihood","hes1M",], breaks = 20, 
     col=rgb(0,0,1,0.3), add=TRUE)

hist(phi2compareMethodsCube["phi2LoocvMse","hes1H",], breaks = 20, 
     col=rgb(1,0,0,0.3))
hist(phi2compareMethodsCube["phi2LoocvLlik","hes1H",], breaks = 20, 
     col=rgb(0,1,0,0.3), add=TRUE)
hist(phi2compareMethodsCube["phi2MarginalLikelihood","hes1H",], breaks = 20, 
     col=rgb(0,0,1,0.3), add=TRUE)

ggobjs <- list()
for(component in c("hes1P", "hes1M", "hes1H")){
  x <- data.frame(log(t(phi2compareMethodsCube[,component,])))
  library(ggplot2);library(reshape2)
  data <- melt(x)
  ggobjs[[component]] <- 
    ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25) + 
    ggplot2::labs(title=component) + xlab("log phi2")
}

ggobjs
ml <- gridExtra::marrangeGrob(ggobjs, nrow=3, ncol=1,
                              top="phi2")
ggsave("phi2-compare-gpFitMethods.pdf", ml, width=8, height = 16)
