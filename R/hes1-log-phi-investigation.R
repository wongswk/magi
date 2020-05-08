library(gpds)
# set up configuration if not already exist ------------------------------------

config <- list(
  seed = 115846061,
  modelName = "Hes1-log"
)

## try phi smoothing ----------------------------

load(paste0("../results/cpp/7param/variablephi-notemper/", config$modelName,"-",config$seed,"-7param-variablephi-notemper.rda"))

xInit <- apply(gpode$xsampled, 2:3, mean)
thetaInit <- colMeans(gpode$theta)
phiNoTemperOptimized <- phiUsed

phiNewInit <- gpds:::solveGpdsRcpp(
  yFull = xInit,
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = xInit,
  thetaInitExogenous = thetaInit,
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = 1,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = 1,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

xdense <- seq(min(xsim$time), max(xsim$time), 0.5)
matplot(x=xsim$time, y=xInit, type="p", pch=1)
matplot(x=xsim$time, y=xsim[,-1], type="p", pch=20, add=TRUE)
gpsmoothFuncList <- list()
for(j in 1:3){
  ynew <- getMeanCurve(xsim.obs$time, xInit[,j], xdense, 
                       t(phiNewInit$phi[,j]), t(pram.true$sigma[j]), kerneltype=config$kernel)
  gpsmoothFuncList[[j]] <- approxfun(xdense, ynew)
  plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time), n=1001,
                lty = 1, col = j, add = TRUE)
}
gpsmoothFuncList <- list()
for(j in 1:3){
  ynew <- getMeanCurve(xsim.obs$time, xInit[,j], xdense, 
                       t(phiExogenous[,j]), t(pram.true$sigma[j]), kerneltype=config$kernel)
  gpsmoothFuncList[[j]] <- approxfun(xdense, ynew)
  plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time), n=1001,
                lty = 2, col = j, add = TRUE)
}


X1 <- na.omit(xsim[,c("time", "X1")])
gpds:::gpsmooth(as.matrix(X1$X1), as.matrix(dist(X1$time)), "generalMatern", sigmaExogenScalar = pram.true$sigma[1], useFrequencyBasedPrior = TRUE)
phiUsed
gpds:::gpsmooth(as.matrix(xInit[,1]), as.matrix(dist(xsim$time)), "generalMatern", sigmaExogenScalar = pram.true$sigma[1], useFrequencyBasedPrior = TRUE)
phiNewInit$phi

yobsThisDim <- X1$X1
priorFactor <- getFrequencyBasedPrior(yobsThisDim)
sigmaExogenScalar = pram.true$sigma[1]
r.nobs <- as.matrix(dist(X1$time))

fn <- function(par) {
  par <- c(par, sigmaExogenScalar)
  marlik <- phisigllikC( par, data.matrix(yobsThisDim), r.nobs, config$kernel)
  penalty <- dnorm(par[2], max(X1$time)*priorFactor["meanFactor"], 
                   max(X1$time)*priorFactor["sdFactor"], log=TRUE)
  -(marlik$value + penalty)
}
gr <- function(par) {
  par <- c(par, sigmaExogenScalar)
  marlik <- phisigllikC( par, data.matrix(yobsThisDim), r.nobs, config$kernel)
  grad <- -as.vector(marlik$grad)
  grad[2] <- grad[2] + (par[2] - max(X1$time)*priorFactor["meanFactor"]) / (max(X1$time)*priorFactor["sdFactor"])^2
  grad[1:2]
}
testthat::expect_equal(gr(c(5,50))[2], (fn(c(5,50+1e-6)) - fn(c(5,50)))/1e-6, tolerance=1e-4)
marlikmap <- optim(rep(1, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2))
marlikmap$par

fn(phiUsed[,1])
fn(phiNewInit$phi[,1])


x1_idx <- which(is.finite(xsim$X1))
gpds:::gpsmooth(as.matrix(xInit[x1_idx,1]), as.matrix(dist(xsim$time[x1_idx])), 
                "generalMatern", sigmaExogenScalar = pram.true$sigma[1], useFrequencyBasedPrior = TRUE)
gpds:::gpsmooth(as.matrix(xInit[x1_idx,1]), as.matrix(dist(xsim$time[x1_idx])), 
                "generalMatern", sigmaExogenScalar = pram.true$sigma[1]/1000, useFrequencyBasedPrior = TRUE)
gpds:::gpsmooth(as.matrix(xInit[x1_idx,1]) + rnorm(length(x1_idx))*pram.true$sigma[1], as.matrix(dist(xsim$time[x1_idx])), 
                "generalMatern", sigmaExogenScalar = pram.true$sigma[1], useFrequencyBasedPrior = TRUE)
yobsThisDim <- xInit[x1_idx,1]
sigmaExogenScalar <- 0.15
# component1 shows that with more smoothed x, bandwidth parameter is big
fn(phiUsed[,1])
fn(phiNewInit$phi[,1])


yobsThisDim <- xInit[,3]
priorFactor <- getFrequencyBasedPrior(yobsThisDim)
sigmaExogenScalar = pram.true$sigma[3]
r.nobs <- as.matrix(dist(xsim$time))
# component3 confirms that bandwidth parameter tends to be too big, gpsmooth is doing its job
fn(phiUsed[,3])
fn(phiNewInit$phi[,3])

# try phi optimization ----------------------------------------------------------------------------------------
samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = xInit,
  thetaInitExogenous = thetaInit,
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = 1/16,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = 1,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

# segfault below
# phi_optimized <- gpds:::optimizePhi(
#   yobsInput = data.matrix(xsim[,-1]), 
#   tvecInput = xsim$time, 
#   fOdeModelInput = hes1logmodel, 
#   sigmaAllDimensionsInput = pram.true$sigma, 
#   priorTemperatureInput = config$priorTemperature, 
#   xInitInput = xInit, 
#   thetaInitInput = thetaInit, 
#   phiInitInput = phiUsed, 
#   missingComponentDim = 2)
