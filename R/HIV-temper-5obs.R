library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 5,  # also set xsim below if changing nobs to 5
    noise = rep(0.3,4),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 400,
    n.iter = 20000,
    burninRatio = 0.5,
    stepSizeFactor = 0.01,
    #filllevel = 3,
    linfillspace = 1,
    t.end = 93,
    modelName = "HIV",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    linearizexInit = TRUE,
    useExoSigma = TRUE,
    useMean = TRUE,
    useBand = TRUE,
    max.epoch = 1
  )
}

config$ndis <- config$t.end / config$linfillspace + 1;

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta = c(0.015, 1.51e-3, 1.11e-3, 4.4e-4, 4.15e-3, 1.1e-3, -2.29e-2, 7.13e-03, 5.68e-04),
  x0 = log(c(3.35e7, 134000, 25000, 9000)),
  sigma=config$noise
)

times <- seq(0, 93, by = 0.25)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::HIVmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = c(70,94,115,139,163) - 70)                # n=5
#xsim <- data.frame(time = c(70,82, 94,106,115,127,139,151,163) - 70) # n=9
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

#xsim <- insertNaN(xsim.obs,config$filllevel)

## Linearly interpolate using fixed interval widths
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}


if (config$useExoSigma) {
  exoSigma = rep(0.001, ncol(xsim)-1) #config$noise
} else {
  exoSigma = numeric(0)
}

if (config$linearizexInit) {  ## linear interpolation for X initializiation
  exoxInit <- sapply(2:ncol(xsim.obs), function(j)
    approx(xsim.obs[, "time"], xsim.obs[, j], xsim[, "time"])$y)
} else {
  exoxInit <- matrix(nrow=0,ncol=0)
}

# cpp inference ----------------------------
HIVmodel <- list(
  name= config$modelName,
  fOde=gpds:::HIVmodelODE,
  fOdeDx=gpds:::HIVmodelDx,
  fOdeDtheta=gpds:::HIVmodelDtheta,
  thetaLowerBound=c(-10,rep(0,5),-10,-10,-10),
  thetaUpperBound=rep(10,9)  
)


## Warm tempering
config$priorTemperatureObs <- 1
config$priorTemperature <- config$ndis / config$nobs  
outDir <- "../results/HIV-temper-5obs/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


### Optimize phi first using equally spaced intervals of 1, i.e., 0,1...,93.
samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = exoxInit,
  odeModel = HIVmodel,
  tvecFull = xsim$time,
  sigmaExogenous = exoSigma,
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = matrix(numeric(0)),
  thetaInitExogenous = matrix(numeric(0)),
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = config$priorTemperatureObs,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = 2,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

phiUsed <- samplesCpp$phi


samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = HIVmodel,
  tvecFull = xsim$time,
  sigmaExogenous = exoSigma,
  phiExogenous = phiUsed,
  xInitExogenous = exoxInit,
  thetaInitExogenous = matrix(numeric(0)),
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = config$priorTemperatureObs,
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

samplesCpp <- samplesCpp$llikxthetasigmaSamples

samplesCpp <- samplesCpp[,,1]

out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim[,-1]))
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(HIVmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))

matplot(xsim$time, xCpp, type="l", add=TRUE)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(HIVmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::HIVmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::HIVmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

#system( paste0("mkdir ", config$seed) )  
save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))
