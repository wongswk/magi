library(magi)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 15,
    noise = rep(0.01, 5), # 0.001 = low noise, 0.01 = high noise
    kernel = "generalMatern",
    loglikflag = "withmeanBand",
    bandsize = 40,
    hmcSteps = 100,
    n.iter = 20001,
    burninRatio = 0.50,
    stepSizeFactor = 0.01,
    linfillspace = 0.5,  # discretization interval width (instead of fill level, since unevenly spaced observations)
    t.end = 100,
    modelName = "PTrans",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    linearizexInit = TRUE,
    useExoSigma = FALSE,
    useMean = TRUE,
    useBand = TRUE,    
    max.epoch = 1
  )
}

# Use this to replicate same seeds used in the paper
args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)[1]
seedlist <- scan("ptrans-noise001-seeds.txt")
config$seed <- seedlist[args]
show(config$seed)


config$ndis <- config$t.end / config$linfillspace + 1;

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,100,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::ptransmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
#matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

## Linearly interpolate using fixed interval widths
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}


if (config$linearizexInit) {
  exoxInit <- sapply(2:ncol(xsim.obs), function(j)
    approx(xsim.obs[, "time"], xsim.obs[, j], xsim[, "time"])$y)
} else {
  exoxInit <- matrix(nrow=0,ncol=0)
}


# cpp inference ----------------------------
ptransmodel <- list(
  name= config$modelName,
  fOde=magi:::ptransmodelODE,
  fOdeDx=magi:::ptransmodelDx,
  fOdeDtheta=magi:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(4,6)
)

outDir <- "../results/ptrans-highnoise/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
config$priorTemperature <- config$ndis / config$nobs
config$priorTemperatureObs <- 1

OursStartTime <- proc.time()[3]

### Optimize phi first using equally spaced intervals of 1, i.e., 0,1...,100.
samplesCpp <- magi:::solveMagiRcpp(
  yFull = exoxInit[xsim$time %in% 0:100,],
  odeModel = ptransmodel,
  tvecFull = 0:100,
  sigmaExogenous = matrix(numeric(0)),
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

samplesCpp <- samplesCpp$llikxthetasigmaSamples
samplesCpp <- samplesCpp[,,1]
out <- samplesCpp[-1,1,drop=FALSE]
sigmaUsed <- tail(out[, 1], ncol(xsim[,-1]))
show(sigmaUsed)

## stabilize phi estimate
samplesCpp <- magi:::solveMagiRcpp(
  yFull = exoxInit[xsim$time %in% 0:100,],
  odeModel = ptransmodel,
  tvecFull = 0:100,
  sigmaExogenous = sigmaUsed,
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

samplesCpp <- samplesCpp$llikxthetasigmaSamples
samplesCpp <- samplesCpp[,,1]
out <- samplesCpp[-1,1,drop=FALSE]
sigmaUsed <- tail(out[, 1], ncol(xsim[,-1]))
show(sigmaUsed)

samplesCpp <- magi:::solveMagiRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = ptransmodel,
  tvecFull = xsim$time,
  sigmaExogenous = sigmaUsed,
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

OursTimeUsed <- proc.time()[3] - OursStartTime

phiUsed <- samplesCpp$phi
samplesCpp <- samplesCpp$llikxthetasigmaSamples

samplesCpp <- samplesCpp[,,1]

out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim[,-1]))
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(ptransmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))

#matplot(xsim$time, xCpp, type="l", add=TRUE)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(ptransmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, magi:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = magi:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, phiUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1],"-fill", config$linfillspace, ".rda"))

