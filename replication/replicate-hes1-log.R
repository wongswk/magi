library(gpds)

outDir <- "../results/hes1log/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) > 0){
  seed <- read.table("R/hes1-log-seeds.txt", header = FALSE)[args, 1]
}else{
  seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9
}

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 33,
    noise = c(0.15,0.15,0.15),
    kernel = "generalMatern",
    seed = seed,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 2e4,
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
# xsim.obs$X1[seq(1,nrow(xsim.obs),by=2)] <- c(0.74, 0.78, 1.86, 1.86, 2.2, 1.93, 1.47, 1.03, 0.36,
#                                              0.88, 1.68, 1.97, 2.15, 1.85, 1.8, 1.47, 0.71)
# xsim.obs$X2[seq(2,nrow(xsim.obs),by=2)] <- c(0.91, 0.82, 0.71, -0.11, 0.08, -0.45, -0.05, 0.2,
#                                              0.88, 1.09, 0.3, 0.35, 0.25, -0.23, -0.51, -0.09)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- insertNaN(xsim.obs,config$filllevel)

# cpp inference ----------------------------
hes1logmodel <- list(
  name="Hes1-log",
  fOde=gpds:::hes1logmodelODE,
  fOdeDx=gpds:::hes1logmodelDx,
  fOdeDtheta=gpds:::hes1logmodelDtheta,
  thetaLowerBound=rep(0,7),
  thetaUpperBound=rep(Inf,7)
)

# heating up with (3, 3, 1) temperature ------------------------------------
#' there are 99 sampled X's (3 x 33) and only 33 actual observations.  
#' Using the logic in our other two models, the default temperatures for Hes1 should be (3, 3, 1).
config$priorTemperature <- 3

samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = matrix(numeric(0)),
  thetaInitExogenous = matrix(numeric(0)),
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = 1,
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


for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(phiUsed[,j], 2), collapse = "; ")
}

save.image(paste0(outDir, config$modelName,"-",config$seed,"-heating.rda"))

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-heating.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)


write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-hes1log_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-hes1log_inferred_theta.csv"))
