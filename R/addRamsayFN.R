# Summarize the results
library(gpds)

config <- list()
config$modelName <- "FN"
config$noise <- rep(0.2, 2)
rdaDir <- "comparison/results/"   ## where ours & Dondel rda saved
outDirWenk <- "comparison/results/"   ## where all the seeds are in separate folders with Wenk's output

seeds <- list.files(outDirWenk, pattern='^\\d+$')  ## get the list of seeds ran

fnmodel <- list(
  fOde=gpds:::fODE,
  fOdeDx=gpds:::fnmodelDx,
  fOdeDtheta=gpds:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)


#### helper function for sample ramsay ------
library(deSolve)
library(CollocInfer)
source("R/helper/utilities.r")
source("R/helper/basic_hmc.R")

sampleRamsay <- function(xsimobs, n.iter, burninRatio){
  FhNdata <- data.matrix(xsimobs[,-1])
  FhNtimes <- xsimobs$time
  
  FhNvarnames <- c("V", "R")
  FhNparnames <- c("a", "b", "c")
  
  colnames(FhNdata) <- FhNvarnames
  
  x0 <- c(-1, 1)
  names(x0) <- FhNvarnames
  
  FhNpars <- c(0.2, 0.2, 3)
  names(FhNpars) <- FhNparnames
  
  fhn.fun <- function(times, x, p, more) {
    dx <- x
    dx[, "V"] <- p["c"] * (x[, "V"] - x[, "V"]^3 / 3 + x[, "R"])
    dx[, "R"] <- -(x[, "V"] - p["a"] + p["b"] * x[, "R"]) / p["c"]
    return(dx)
  }
  
  FhNn <- length(FhNtimes)
  
  # set up B-splines.
  FhNrange <- c(0, 20)
  breaks <- seq(0, 20, 0.1)
  FhNbasis <- create.bspline.basis(range = FhNrange, norder = 4, breaks = breaks)
  
  FhNfdPar <- fdPar(FhNbasis, int2Lfd(2), 1)
  lambda <- 1000  # start pars at truth so this OK, Ramsey method's estimates converge as lambda --> infty
  # initial smooth
  DEfd0 <- smooth.basis(FhNtimes, FhNdata, FhNfdPar)$fd
  coefs0 <- DEfd0$coef
  
  profile.obj <- LS.setup(pars = FhNpars, fn = fhn.fun, lambda = lambda,
                          times = FhNtimes, coefs = coefs0, basisvals = FhNbasis)
  
  proc <- profile.obj$proc
  lik <- profile.obj$lik
  
  Ores2.2 <- Profile.LS(fhn.fun, FhNdata, FhNtimes, FhNpars, coefs0, FhNbasis, lambda)    # Ramsey optimization, use this as starting values for sampling
  
  Ores2.2$pars
  
  ################################
  
  
  # assume current pars (curpars) defined
  Cblock <- function(q, grad=FALSE) {
    
    data<- FhNdata
    times <- FhNtimes
    pars <- curpars
    
    coefs2 = matrix(q, ncol(lik$bvals), length(q)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2, proc$bvals, pars, proc$more)
    
    ret <- -f
    
    if(grad) {
      g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals, 
                                              pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars, proc$more)
      g = as.vector(g)    
      attr(ret,"grad") <- -g
    }
    
    return(ret)
  }
  
  
  # assume current spline coefficients (curcoefs) are defined
  parblock <- function(q, grad=FALSE) {
    
    data<- FhNdata
    times <- FhNtimes
    coefs <- curcoefs
    pars <- q
    
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2, proc$bvals, pars, proc$more)
    
    ret <- -f
    
    if(grad) {
      attr(ret,"grad") <- -SplineCoefsDP(coefs, times, data, lik, proc, pars)
    }
    
    return(ret)
  }
  
  
  #### actual sampling
  
  th.pars <- matrix(NA,n.iter,3)  
  th.C <- matrix(NA,n.iter,406) 
  th.pars[1,] <- Ores2.2$pars
  th.C[1,] <- c(Ores2.2$coefs[,1],Ores2.2$coefs[,2])
  curpars <- th.pars[1,]
  curcoefs <- th.C[1,]
  lliklist <- c()
  lliklist[1] <- parblock(Ores2.2$pars)
  accepts <- matrix(NA,n.iter,2)
  accepts[1,] <- c(1,1)
  
  
  for (t in 2:n.iter) {
    
    if (t %% 100 == 0) { show(c(t, lliklist[t-1], apply(accepts[1:(t-1),],2,mean))) }
    
    # Update pars
    foo <- basic_hmc(parblock, step=runif(1,0.015,0.03)* c(1,1,2.5), nsteps= 20, initial=th.pars[t-1,], return.traj = T)
    th.pars[t,] <- foo$final
    curpars <- th.pars[t,]
    accepts[t,1] <- foo$acc
    
    # Update C
    foo <- basic_hmc(Cblock, step=runif(1,0.005,0.01), nsteps= 20, initial=th.C[t-1,], return.traj = T)
    th.C[t,] <- foo$final
    curcoefs <- th.C[t,]
    accepts[t,2] <- foo$acc
    
    lliklist[t] <- foo$lpr
    
  }
  
  burnin <- as.integer(n.iter * burninRatio)
  
  gpode <- list(theta=th.pars[-(1:burnin), ],
                xsampled=array(th.C[-(1:burnin),],
                               dim=c(n.iter-burnin, ncol(th.C) / ncol(FhNdata), ncol(FhNdata))),
                lglik=lliklist[-(1:burnin)])
  
  return(gpode)
}


#### Grab the data for each seed and save in a list
args <- commandArgs(trailingOnly = TRUE)
j <- as.numeric(args)
i <- seeds[j]

load(paste0(rdaDir, config$modelName,"-",i,"-noise", config$noise[1], ".rda"))
configRamsay <- list(
  n.iter=20000,
  burninRatio=0.4,
  modelName="FN",
  seed=config$seed,
  noise=config$noise,
  nobs=config$nobs
)
gpodeRamsay <- sampleRamsay(na.omit(xsim), configRamsay$n.iter, configRamsay$burninRatio)
saveRDS(gpodeRamsay, paste0(rdaDir, configRamsay$modelName,"-",i,"-noise", configRamsay$noise[1], "-gpodeRamsay.rds"))

gpodeRamsay$sigma <- matrix(pram.true$sigma, nrow=length(gpodeRamsay$lglik), ncol=length(pram.true$sigma), byrow = TRUE)
xsimRamsay <- data.frame(time=seq(0, 20, 0.1), V=NA, R=NA)
xsimRamsay$V[match(na.omit(xsim)$time, xsimRamsay$time)] <- na.omit(xsim)[,2]
xsimRamsay$R[match(na.omit(xsim)$time, xsimRamsay$time)] <- na.omit(xsim)[,3]
gpodeRamsay$xsampled <- gpodeRamsay$xsampled[,2:202,]
gpodeRamsay$fode <- sapply(1:length(gpodeRamsay$lglik), function(t) 
  with(gpodeRamsay, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpodeRamsay$fode <- aperm(gpodeRamsay$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0(rdaDir, configRamsay$modelName,"-",configRamsay$seed,"-noise", configRamsay$noise[1], "-gpodeRamsay.pdf"), 
  xtrue, dotxtrue, xsimRamsay, gpodeRamsay, pram.true, configRamsay, odemodel)
