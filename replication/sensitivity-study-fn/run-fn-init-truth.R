library(magi)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) > 0){
  filllevel <- args[1]
  seed <- scan("fn-seeds.txt")[args[2]]
  nobs_keep <- args[3]
}else{
  seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9
  filllevel <- 2
  nobs_keep <- 41
}

flagInitX <- "truth"

outDir <- paste0("../results/fn-init-truth-fill", filllevel, "-nobs", nobs_keep, "/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = nobs_keep,
    noise = c(0.2, 0.2),
    kernel = "generalMatern",
    seed = seed,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = filllevel,
    t.end = 20,
    modelName = "FN",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    max.epoch = 1
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperatureObs <- 1

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

times <- seq(0,config$t.end,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::fnmodelODE(parameters, t(state))))
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

fnmodel <- list(
  fOde=magi:::fODE,
  fOdeDx=magi:::fnmodelDx,
  fOdeDtheta=magi:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

config$priorTemperature <- config$ndis / config$nobs  
config$priorTemperatureObs <- 1

if(flagInitX == "linear"){
  xInitExogenous <- data.matrix(xsim[,-1])
  for (j in 1:(ncol(xsim)-1)){
    xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
  }
  thetaInitExogenous <- matrix(nrow=0,ncol=0)
  sigmaExogenous <- numeric(0)
}else if(flagInitX == "truth"){
  xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
  thetaInitExogenous <- pram.true$theta
  sigmaExogenous <- config$noise
}

OursStartTime <- proc.time()[3]

samplesCpp <- magi:::solveMagiRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = fnmodel,
  tvecFull = xsim$time,
  sigmaExogenous = sigmaExogenous,
  phiExogenous = matrix(nrow=0,ncol=0),
  xInitExogenous = xInitExogenous,
  thetaInitExogenous = thetaInitExogenous,
  muExogenous = matrix(nrow=0,ncol=0),
  dotmuExogenous = matrix(nrow=0,ncol=0),
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
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(fnmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))

matplot(xsim$time, xCpp, type="l", add=TRUE)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(fnmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, magi:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = magi:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(phiUsed[,j], 2), collapse = "; ")
}

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_theta.csv"))
