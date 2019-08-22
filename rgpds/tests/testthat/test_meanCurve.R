#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = c(0.15, 0.07) * 2,
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 1e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = 2,
    modelName = "FN",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = FALSE,
    useGPmean = TRUE,
    forseTrueMean = FALSE,
    phase2 = FALSE,
    temperPrior = TRUE,
    phase3 = FALSE,
    max.epoch = 10,
    epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1]
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
if(config$temperPrior){
  config$priorTemperature <- config$ndis / config$nobs  
}else{
  config$priorTemperature <- 1
}

# initialize global parameters, true x, simulated x ----------------------------
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
  config$seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9 # random seed on cluster
}else{
  baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                              "-ndis",ndis,"-",ifelse(config$temperPrior, "temperPrior", "unitHeatPrior"),"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,20,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,20,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

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

# GPsmoothing: marllik+fftNormalprior for phi-sigma ----------------------------
cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
  
  desiredMode <- priorFactor["meanFactor"]
  betaRate <- uniroot(function(betaRate) pgamma(1, 1 + desiredMode*betaRate, betaRate)-0.95,
                      c(1e-3, 1e3))$root
  alphaRate <- 1 + desiredMode*betaRate
  
  fn <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                     max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
    # penalty <- dgamma(par[2], alphaRate, betaRate/max(xsim.obs$time), log=TRUE)
    # penalty <- 0
    -(marlik$value + penalty)
  }
  gr <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad)
    penalty <- (par[2] - max(xsim.obs$time)*priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
    # penalty <- ((alphaRate-1)/par[2] - betaRate/max(xsim.obs$time))
    # penalty <- 0
    grad[2] <- grad[2] + penalty
    grad
  }
  testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-3)
  marlikmap <- optim(c(sd(xsim.obs[,1+j])/2, max(xsim.obs$time)/2, sd(xsim.obs[,1+j])/2), 
                     fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, Inf, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}

cursigma
curphi

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})


for(j in 1:(ncol(xsim)-1)){
  ydyR <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                       t(curphi[,j]), t(cursigma[j]), 
                       kerneltype=config$kernel, deriv = TRUE)
  ydyC <- gpds:::calcMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time,
                               curphi[,j,drop=FALSE], cursigma[j],
                               "generalMatern", TRUE)
  testthat::expect_equal(ydyR[[1]], ydyC[,,1], check.attributes = FALSE)
}


gpsmoothFuncList <- list()
for(j in 1:(ncol(xsim)-1)){
  ynew <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                       t(curphi[,j]), t(cursigma[j]), kerneltype=config$kernel)
  gpsmoothFuncList[[j]] <- approxfun(xsim$time, ynew)
  plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time),
                lty = 2, col = j, add = TRUE)
}
cursigma
# MCMC starting value ----------------------------------------------------------
yobs <- data.matrix(xsim[,-1])

xsimInit <- xsim
for(j in 1:(ncol(xsim)-1)){
  nanId <- which(is.na(xsimInit[,j+1]))
  xsimInit[nanId,j+1] <- gpsmoothFuncList[[j]](xsimInit$time[nanId])
}
matplot(xsimInit$time, xsimInit[,-1], type="p", pch=2, add=TRUE)
xInit <- data.matrix(xsimInit[,-1])

thetaInit <- rep(1, length(pram.true$theta))

thetaoptim <- function(xInit, thetaInit, curphi, cursigma){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    -xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "FN" )$value
  }
  gr <- function(par) {
    -as.vector(xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "FN" )$grad[-(1:length(xInit))])
  }
  marlikmap <- optim(c(thetaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  thetaInit[] <- marlikmap$par
  list(thetaInit = thetaInit)
}

thetamle <- thetaoptim(xInit, thetaInit, curphi, cursigma)
