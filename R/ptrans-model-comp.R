library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 15,
    noise = rep(0.01, 5), # 0.001 = low noise, 0.01 = high noise
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 40,
    hmcSteps = 100,
    n.iter = 20001,
    n.iter.Wenk = 300000,
    n.iter.Dondel = 300000,
    burninRatio = 0.50,
    stepSizeFactor = 0.01,
    #filllevel = 3,
    linfillspace = 0.1,
    t.end = 100,
    modelName = "PTrans",
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

#config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$ndis <- config$t.end / config$linfillspace + 1;
if(config$temperPrior){
  config$priorTemperature <- config$ndis / config$nobs  
}else{
  config$priorTemperature <- 1
}

# initialize global parameters, true x, simulated x ----------------------------
outDir <- "../results/cpp/"
system(paste("mkdir -p", outDir))
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,100,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::ptransmodelODE(parameters, t(state))))
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

#xsim <- insertNaN(xsim.obs,config$filllevel)
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}

if (config$useExoSigma) {
  exoSigma = rep(0.001, ncol(xsim)-1)
} else {
  exoSigma = numeric(0)
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
  fOde=gpds:::ptransmodelODE,
  fOdeDx=gpds:::ptransmodelDx,
  fOdeDtheta=gpds:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(4,6)
)

OursStartTime <- proc.time()[3]

samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = ptransmodel,
  tvecFull = xsim$time,
  sigmaExogenous = exoSigma,
  phiExogenous = matrix(numeric(0)),
  xInitExogenous = exoxInit,
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

OursTimeUsed <- proc.time()[3] - OursStartTime

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
  with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

### Now run Dondel
library(deGradInfer)

VG_func <- function(t, X, params) {
  MM <- ((params[5]*X[,5])/(params[6]+X[,5]))# Michaelis-Menten term
  dxdt <- cbind( - params[1]*X[,1] - params[2]*X[,1]*X[,3] +params[3]*X[,4],# S
                 params[1]*X[,1],# dS
                 - params[2]*X[,1]*X[,3] + params[3]*X[,4] +MM,# R
                 params[2]*X[,1]*X[,3] - params[3]*X[,4] -params[4]*X[,4],# RS
                 params[4]*X[,4] - MM# Rpp
  )
  return(dxdt)
}

# Copy data
timeTest <- xsim.obs[,1]  #c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100)
dataTest <- as.matrix(xsim.obs[,2:ncol(xsim)])

# define our own prior with parameter bounds
log_prior <- function(params) {
  return(c(dunif(params[1],0,4,log=TRUE),
           dunif(params[2],0,4,log=TRUE),
           dunif(params[3],0,4,log=TRUE),
           dunif(params[4],0,4,log=TRUE),
           dunif(params[5],0,4,log=TRUE),
           dunif(params[6],0,4,log=TRUE)))
}

DondelStartTime <- proc.time()[3]

agm.result =agm(data=dataTest,time=timeTest,ode.system=VG_func, numberOfParameters=length(pram.true$theta),
                noise.sd = config$noise[1], logPrior = log_prior, maxIterations = config$n.iter.Dondel, showProgress = TRUE)

DondelTimeUsed <- proc.time()[3] - DondelStartTime

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
  with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)

save(agm.result, file=paste0(outDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".rda"))

### Wenk
#### Set up a place to store Python Wenk temporary output  ### Remember to set Noise and Iterations in Python!

WenkStartTime <- proc.time()[3]

system( paste0("mkdir ", config$seed) )  
Sys.setenv(PYTHONPATH="../comparison/")
system( paste0("cd ", config$seed, "; python3 ../comparison/FGPGM/mainFiles/ProteinTransduction/createExperiments.py" ))
## copy our xsim.obs to observations.csv so same dataset is used
write.table(as.matrix(xsim.obs[,2:ncol(xsim)]), row.names = F, col.names = F, file=paste0(config$seed, "/observations.csv"))

system( paste0("cd ", config$seed, "; python3 ../comparison/FGPGM/mainFiles/ProteinTransduction/getHyperparams.py" ))
system( paste0("cd ", config$seed, "; python3 ../comparison/FGPGM/mainFiles/ProteinTransduction/doFGPGM.py" ))

WenkTimeUsed <- proc.time()[3] - WenkStartTime

cat(OursTimeUsed, DondelTimeUsed, WenkTimeUsed, file=paste0(outDir, config$modelName, "-time-", config$seed,"-noise", config$noise[1], ".txt"))

dataDir <- paste0(config$seed, "/")

## use our plot codes
samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))
lliklist <- as.vector(as.matrix(read.table(paste0(dataDir, "lliklist.csv") )))

llikId <- 0  ### llik is in its own file
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim.obs[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(ptransmodel$thetaLowerBound))
#sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))  ## no sigma sampled

## Un-standardize their X's and untransform thetas
samplesCpp[,xId] <- t(apply(samplesCpp[,xId],1, function(x) xSD*x + xMean))
samplesCpp[,thetaId] <- t(apply(samplesCpp[,thetaId],1, function(x) x * 10^thetaMag))

burnin <- as.integer(config$n.iter.Wenk*config$burninRatio)
gpode <- list(theta= samplesCpp[-(1:burnin), thetaId],
              xsampled=array(samplesCpp[-(1:burnin), xId],
                             dim=c(nrow(samplesCpp)-burnin, nrow(xsim.obs), ncol(xsim)-1)),
              lglik=  lliklist[-(1:burnin)], ###samplesCpp[llikId,-(1:burnin)],
              sigma=  matrix(0, nrow = nrow(samplesCpp) - burnin, ncol = ncol(xsim)-1)) # t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-Wenk-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)





