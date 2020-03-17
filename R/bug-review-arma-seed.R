#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 33,
    noise = c(0.15,0.15,0.15),
    kernel = "generalMatern",
    seed = 480330981, # seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 5e2+2,
    burninRatio = 0.50,
    stepSizeFactor = 0.001,
    filllevel = 0,
    modelName = "Hes1-log",
    async = TRUE,
    max.epoch = 1,
    useMean = TRUE,
    useBand = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = TRUE
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs

# initialize global parameters, true x, simulated x ----------------------------
outDir <- "../results/cpp/"
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = log(c(1.438575, 2.037488, 17.90385)),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1logmodelODE(parameters, t(state))))
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
xsim$X3 <- NaN
xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
if(config$async){
  xsim.obs$X1[seq(2,nrow(xsim.obs),by=2)] <- NaN
  xsim.obs$X2[seq(1,nrow(xsim.obs),by=2)] <- NaN
}
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- insertNaN(xsim.obs,config$filllevel)

# cpp inference ----------------------------
hes1logmodelODE_restricted <- function(theta, x){
  gpds:::hes1logmodelODE(c(theta[1:5], 20, theta[6]), x)
}

hes1logmodelDx_restricted <- function(theta, x){
  gpds:::hes1logmodelDx(c(theta[1:5], 20, theta[6]), x)
}

hes1logmodelDtheta_restricted <- function(theta, x){
  gpds:::hes1logmodelDtheta(c(theta[1:5], 20, theta[6]), x)[,c(1:5,7),]
}

hes1logmodel <- list(
  # name="Hes1-log",
  fOde=hes1logmodelODE_restricted,
  fOdeDx=hes1logmodelDx_restricted,
  fOdeDtheta=hes1logmodelDtheta_restricted,
  thetaLowerBound=c(rep(0,6)),
  thetaUpperBound=c(rep(Inf,6))
)

# use the phi from hes1-log-async-partialobs-cpp.R with seed 1365546660
phiExogenous <- rbind(
  c(2.1282, 0.3872, 0.4645),
  c(65.1345, 35.9814, 22.4503)
)

hes1logmodelcpp <- list(
  name="Hes1-log-fixf-noname",
  fOde=gpds:::hes1logmodelODEfixf,
  fOdeDx=gpds:::hes1logmodelDxfixf,
  fOdeDtheta=gpds:::hes1logmodelDthetafixf,
  #thetaLowerBound=rep(0,7),
  #thetaUpperBound=rep(Inf,7)
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(Inf,6)
)

param_restricted <- pram.true
param_restricted$theta <- c(param_restricted$theta[1:5], 0.3)

test_theta <- param_restricted$theta * exp(rnorm(6)/10)
noise <- rnorm(length(data.matrix(xtrue)[,-1]))*4
out_r <- hes1logmodel$fOde(test_theta, data.matrix(xtrue)[,-1]+noise)
out_cpp <- hes1logmodelcpp$fOde(test_theta, data.matrix(xtrue)[,-1]+noise)
testthat::expect_equal(out_r, out_cpp)

out_r <- hes1logmodel$fOdeDx(test_theta, data.matrix(xtrue)[,-1]+noise)
out_cpp <- hes1logmodelcpp$fOdeDx(test_theta, data.matrix(xtrue)[,-1]+noise)
testthat::expect_equal(out_r, out_cpp)

out_r <- hes1logmodel$fOdeDtheta(test_theta, data.matrix(xtrue)[,-1]+noise)
out_cpp <- hes1logmodelcpp$fOdeDtheta(test_theta, data.matrix(xtrue)[,-1]+noise)
testthat::expect_equal(out_r, out_cpp)

set.seed(config$seed)
runif(2)

# cpp has bug
samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodelcpp,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = phiExogenous,
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

# r is fine
samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = phiExogenous,
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
