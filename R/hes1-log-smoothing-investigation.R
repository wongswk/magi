library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 33,
    noise = c(0.15,0.15,0.1),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 2e2,
    burninRatio = 0.50,
    stepSizeFactor = 0.01,
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

hes1logmodel <- list(
  name="Hes1-log",
  fOde=gpds:::hes1logmodelODE,
  fOdeDx=gpds:::hes1logmodelDx,
  fOdeDtheta=gpds:::hes1logmodelDtheta,
  thetaLowerBound=rep(0,7),
  thetaUpperBound=rep(Inf,7)
)

## tempered with warm start and update phi from warm start ----------------------------

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

samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = phiNewInit$phi,
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

phiUsed <- samplesCpp$phi
samplesCpp <- samplesCpp$llikxthetasigmaSamples

samplesCpp <- samplesCpp[,,1]
out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim[,-1]))
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(hes1logmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))

matplot(xsim$time, xCpp, type="l", add=TRUE, lty=2)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(hes1logmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1logmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1logmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- "../results/cpp/7param/variablephi-temper-warmstart-updatephi/"
system(paste("mkdir -p", outDir))

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(phiUsed[,j], 2), collapse = "; ")
}

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-7param-variablephi-temper-warmstart-updatephi.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save.image(paste0(outDir, config$modelName,"-",config$seed,"-7param-variablephi-temper-warmstart-updatephi.rda"))
