#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = c(0.15, 0.07) * 2,
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 101,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = 2,
    modelName = "FN",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    max.epoch = 1
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
if(config$temperPrior){
  config$priorTemperature <- config$ndis / config$nobs
}else{
  config$priorTemperature <- 1
}

if(config$loglikflag == "withmeanBand"){
  config$useMean = TRUE
  config$useBand = TRUE
}else if(config$loglikflag == "band"){
  config$useMean = FALSE
  config$useBand = TRUE
}else if(config$loglikflag == "withmean"){
  config$useMean = TRUE
  config$useBand = FALSE
}else if(config$loglikflag == "usual"){
  config$useMean = FALSE
  config$useBand = FALSE
}

# initialize global parameters, true x, simulated x ----------------------------
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


# cpp inference ----------------------------
fnmodel <- list(
  fOde=gpds:::fODE,
  fOdeDx=gpds:::fnmodelDx,
  fOdeDtheta=gpds:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf)
)

samplesCpp <- gpds:::solveGpds(
  yFull = data.matrix(xsim[,-1]),
  odeModel = fnmodel,
  tvecFull = xsim$time,
  sigmaExogenous = numeric(0),
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = matrix(numeric(0)),
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = config$n.iter,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=2)
stopifnot(abs(sum(out[,1])*1e5 - 6879957.07974693) < 1)
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + 3), 1]
sigmaCpp <- tail(out[, 1], 2)

matplot(xsim$time, xCpp, type="l", add=TRUE)

phiExogenous = cbind(c(2.24, 1.64), c(0.65, 2.93))

samplesCpp <- gpds:::solveGpds(
  yFull = data.matrix(xsim[,-1]),
  odeModel = fnmodel,
  tvecFull = xsim$time,
  sigmaExogenous = c(0.30, 0.15),
  xInitExogenous = xCpp,
  muExogenous = matrix(0, nrow=nrow(xsim), ncol=2),
  dotmuExogenous = matrix(0, nrow=nrow(xsim), ncol=2),
  phiExogenous = phiExogenous,
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = config$n.iter,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

samplesCpp <- gpds:::solveGpds(
  yFull = data.matrix(xsim[,-1]),
  odeModel = fnmodel,
  tvecFull = xsim$time,
  sigmaExogenous = numeric(0),
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = matrix(numeric(0)),
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = config$n.iter,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = 4,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)


# R inference ----------------------------

## GPsmoothing: marllik+fftNormalprior for phi-sigma; CHECKS OUT
tvec.nobs <- xsim.obs$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j]) # here has discrepancy because quantile function in R and c++ are different
  
  desiredMode <- priorFactor["meanFactor"]
  
  fn <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                     max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
    -(marlik$value + penalty)
  }
  gr <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad)
    penalty <- (par[2] - max(xsim.obs$time)*priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
    grad[2] <- grad[2] + penalty
    grad
  }
  testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-3)
  if(j == 1){
    phisigCpp <- c(2.2433, 1.6356, 0.3352574)
  }else if(j == 2){
    phisigCpp <- c(0.6527, 2.9276, 0.1472876)
  }
  marlikmap <- optim(c(phisigCpp), 
                     fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, Inf, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
  testthat::expect_equal(marlikmap$par, phisigCpp, tolerance=1e-5)
}

cursigma
curphi

tvec.full <- xsim$time
foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

if(config$useMean){
  for(j in 1:(ncol(xsim)-1)){
    ydy <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                        t(curphi[,j]), t(cursigma[j]), 
                        kerneltype=config$kernel, deriv = TRUE)
    curCov[[j]]$mu <- as.vector(ydy[[1]])
    curCov[[j]]$dotmu <- as.vector(ydy[[2]])
  }
}

## initXmudotmu
# can explicitly export the cov to see further
for(j in 1:(ncol(xsim)-1)){
  testthat::expect_equal(curCov[[j]]$mu, xCpp[,j], tolerance=1e-5)
}
xInit <- cbind(curCov[[1]]$mu, curCov[[2]]$mu)

## initTheta
thetaoptim <- function(xInit, thetaInit, cursigma){
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
thetaInit <- rep(1, 3)
yobs <- data.matrix(xsim[,-1])
thetamle <- thetaoptim(xInit, thetaInit, cursigma)
testthat::expect_equal(thetamle$thetaInit, thetaCpp, tolerance=5e-5)
thetaInit <- thetamle$thetaInit

## hmc sampler
sigmaInit <- cursigma
xthetasigmaInit <- c(xInit, thetaInit, sigmaInit)
stepLowInit <- rep(config$stepSizeFactor/config$hmcSteps, length(xthetasigmaInit))

xId <- 1:length(xInit)
thetaId <- (max(xId)+1):(max(xId)+length(thetaInit))
sigmaId <- (max(thetaId)+1):(max(thetaId)+length(sigmaInit))

xthetasigamSingleSampler <- function(xthetasigma, stepSize) 
  xthetasigmaSample(yobs, curCov, xthetasigma[sigmaId], xthetasigma[c(xId, thetaId)], 
                    stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
                    priorTemperature = config$priorTemperature, modelName = config$modelName)

chainSamplesOut <- chainSampler(config, xthetasigmaInit, xthetasigamSingleSampler, stepLowInit, verbose=TRUE)
