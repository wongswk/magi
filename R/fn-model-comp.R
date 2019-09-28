#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
#where to save results
PROJECT_DIR = getwd()
outDir <- "comparison/results/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = c(0.2, 0.2),
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    n.iter.Wenk = 300000,
    n.iter.Dondel = 300000,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = 2,
    t.end = 20,
    modelName = "FN",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    useExoSigma = TRUE,
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

times <- seq(0,config$t.end,length=241)

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

if (config$useExoSigma) {
  exoSigma = config$noise
} else {
  exoSigma = numeric(0)
}

# cpp inference ----------------------------
fnmodel <- list(
  fOde=gpds:::fODE,
  fOdeDx=gpds:::fnmodelDx,
  fOdeDtheta=gpds:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = fnmodel,
  tvecFull = xsim$time,
  sigmaExogenous = exoSigma,
  phiExogenous = matrix(nrow=0,ncol=0),
  xInitExogenous = matrix(nrow=0,ncol=0),
  muExogenous = matrix(nrow=0,ncol=0),
  dotmuExogenous = matrix(nrow=0,ncol=0),
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
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

### Now run Dondel
library(deGradInfer)

FN_func <- function(t, x, params) {
  gpds:::fnmodelODE(params, x)
}


# Copy data
timeTest <- xsim.obs[,1]
dataTest <- as.matrix(xsim.obs[,2:ncol(xsim)])

# define our own prior with parameter bounds
log_prior <- function(params) {
  return(c(dunif(params[1],0,10,log=TRUE),
           dunif(params[2],0,10,log=TRUE),
           dunif(params[3],0,10,log=TRUE)))
}

agm.result =agm(data=dataTest,time=timeTest,ode.system=FN_func, numberOfParameters=length(pram.true$theta),
                noise.sd = config$noise[1], logPrior = log_prior, maxIterations = config$n.iter.Dondel, showProgress = TRUE)

#### Use our plotting codes
config$n.iter.Dondel <- config$n.iter.Dondel / 25
x_means <- apply(xsim.obs[,-1],2,mean)

xsampled <- (sapply(agm.result$x.samples, function(x) x, simplify="array"))
xsampled <- aperm(apply(xsampled, c(1,2), function(x) x + x_means), c(2,3,1))
thetasampled <- agm.result$posterior.samples
lglik <- agm.result$ll

burnin <- as.integer(config$n.iter.Dondel*config$burninRatio)

gpode <- list(theta= thetasampled[-(1:burnin),],
              xsampled= xsampled[-(1:burnin),,],
              lglik=  lglik[-(1:burnin)],
              sigma=  matrix(0, nrow = config$n.iter-burnin, ncol = ncol(xsim)-1)) # not sampled in this method
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)

save(agm.result, file=paste0(outDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".rda"))

### Wenk
#### Set up a place to store Python Wenk temporary output  ### Remember to set Noise and Iterations in Python!
setwd(outDir)
system( paste0("mkdir ", config$seed) )  
Sys.setenv(PYTHONPATH=paste0(PROJECT_DIR, "/comparison/"))
system( paste0("cd ", config$seed, "; python3 ", PROJECT_DIR, "/comparison/FGPGM/mainFiles/FitzHughNagumo/createExperiments.py", 
               " --tEnd ", config$t.end, " --dt ", config$t.end / (config$nobs - 1), " --obsNoiseStd ", config$noise[1]))
## copy our xsim.obs to observations.csv so same dataset is used
write.table(as.matrix(xsim.obs[,2:ncol(xsim)]), row.names = F, col.names = F, file=paste0(config$seed, "/observations.csv"))

system( paste0("cd ", config$seed, "; python3 ", PROJECT_DIR, "/comparison/FGPGM/mainFiles/FitzHughNagumo/getHyperparams.py" ))
system( paste0("cd ", config$seed, "; python3 ", PROJECT_DIR, "/comparison/FGPGM/mainFiles/FitzHughNagumo/doFGPGM.py",
               " --nSamples ", sprintf("%d", config$n.iter.Wenk)))

dataDir <- paste0(config$seed, "/")

## use our plot codes
samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))
lliklist <- as.vector(as.matrix(read.table(paste0(dataDir, "lliklist.csv") )))
sigmaWenk <- list()
for(j in 1:(ncol(xsim)-1)){
  sigmaWenk[[j]] <- read.table(paste0(dataDir, "/hyperparams/sigma", j - 1, ".csv") )
}
sigmaWenk <- unlist(sigmaWenk)

llikId <- 0  ### llik is in its own file
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim.obs[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(fnmodel$thetaLowerBound))
#sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))  ## no sigma sampled

## Un-standardize their X's and untransform thetas
samplesCpp[,xId] <- t(apply(samplesCpp[,xId],1, function(x) xSD*x + xMean))
samplesCpp[,thetaId] <- t(apply(samplesCpp[,thetaId],1, function(x) x * 10^thetaMag))

burnin <- as.integer(config$n.iter.Wenk*config$burninRatio)
gpode <- list(theta= samplesCpp[-(1:burnin), thetaId],
              xsampled=array(samplesCpp[-(1:burnin), xId],
                             dim=c(nrow(samplesCpp)-burnin, nrow(xsim.obs), ncol(xsim)-1)),
              lglik=  lliklist[-(1:burnin)], ###samplesCpp[llikId,-(1:burnin)],
              sigma=  matrix(sigmaWenk, nrow = nrow(samplesCpp) - burnin, ncol = ncol(xsim)-1, byrow = TRUE)) # t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

setwd(PROJECT_DIR)
gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-Wenk-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
