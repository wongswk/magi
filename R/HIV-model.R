#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
source("R/HIV-helper-functions.R")
# set up configuration if not already exist ------------------------------------
#if(!exists("config")){
  config <- list(
    nobs = 11,
    noise = rep(0.5,4),
    kernel = "generalMatern",
    seed = 396033147,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 200,
    n.iter = 40000,
    burninRatio = 0.5,
    stepSizeFactor = 0.1,
    filllevel = 1,
    modelName = "HIV",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = FALSE,
    useGPmean = TRUE
  )
#}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
# config$priorTemperature[2] <- 1e12

# initialize global parameters, true x, simulated x ----------------------------
if(grepl("/ufrc/",getwd())){
  baseDir <- "/ufrc/swong/swkwong/gpds/results/" # tmp folder on cluster 
}else if(grepl("/shihaoyang/",getwd())){
  baseDir <- "~/Workspace/DynamicSys/results/"  
}else{
  baseDir <- "c:/data/projects/DynamicSys/results/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                              "-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.015, 1.51e-3, 1.11e-3, 4.4e-4, 4.15e-3, 1.1e-3, -2.29e-2, 7.13e-03, 5.68e-04),
  xinit = log(c(3.2e7, 134000, 26000, 10000)),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 100, by = 0.25)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::HIVmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$xinit, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- xtrue

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- insertNaN(xsim.obs,config$filllevel)

tvec.full <- xsim$time
tvec.nobs <- xsim.obs$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

# GPsmoothing: marllik+fftNormalprior for phi-sigma ----------------------------
eval(phiAllMethodsExpr)
cursigma <- cursigmaMarllikFftprior
curphi <- curphiMarllikFftprior

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

if(config$useGPmean){
  for(j in 1:(ncol(xsim)-1)){
    ydy <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                        t(curphi[,j]), t(cursigma[j]), 
                        kerneltype=config$kernel, deriv = TRUE)
    curCov[[j]]$mu <- as.vector(ydy[[1]])
    curCov[[j]]$dotmu <- as.vector(ydy[[2]])
  }
}

gpsmoothFuncList <- list()
for(j in 1:4){
  ynew <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                       t(curphi[,j]), t(cursigma[j]), kerneltype=config$kernel)
  gpsmoothFuncList[[j]] <- approxfun(xsim$time, ynew)
  plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time),
                lty = 2, col = j, add = TRUE)
}
cursigma
# MCMC starting value ----------------------------------------------------------
yobs <- data.matrix(xsim[,-1])


xsimInit <- xsim
for(j in 1:4){
  nanId <- which(is.na(xsimInit[,j+1]))
  xsimInit[nanId,j+1] <- gpsmoothFuncList[[j]](xsimInit$time[nanId])
}
matplot(xsimInit$time, xsimInit[,-1], type="p", pch=2, add=TRUE)
xInit <- data.matrix(xsimInit[,-1])

thetaInit <- rep(1, length(pram.true$theta))
thetamle <- thetaoptim(xInit, thetaInit, curphi, cursigma)
thetaInit <- thetamle$thetaInit
sigmaInit <- cursigma

if(config$startXAtTruth){
  xInit <- sapply(xtrueFunc, function(f) f(xsim$time))  
}
if(config$startThetaAtTruth){
  thetaInit <- pram.true$theta  
}
if(config$startSigmaAtTruth){
  sigmaInit <- pram.true$sigma  
}

xthetasigmaInit <- c(xInit, thetaInit, sigmaInit)
stepLowInit <- rep(0.000035, length(xthetasigmaInit))
stepLowInit <- stepLowInit*config$stepSizeFactor


# HMC sampler for x, theta, sigma ----------------------------------------------
xId <- 1:length(xInit)
thetaId <- (max(xId)+1):(max(xId)+length(thetaInit))
sigmaId <- (max(thetaId)+1):(max(thetaId)+length(sigmaInit))

xthetasigamSingleSampler <- function(xthetasigma, stepSize) 
  xthetasigmaSample(yobs, curCov, xthetasigma[sigmaId], xthetasigma[c(xId, thetaId)], 
                    stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
                    priorTemperature = config$priorTemperature, modelName = config$modelName)
# stepLowInit[sigmaId] <- 0
chainSamplesOut <- chainSampler(config, xthetasigmaInit, xthetasigamSingleSampler, stepLowInit, verbose=TRUE)

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), thetaId],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), xId], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = chainSamplesOut$xth[-(1:burnin), sigmaId, drop=FALSE],
              phi = curphi)
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::HIVmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::HIVmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
names(philist) <- paste0("phi", 1:length(philist))
configWithPhiSig <- c(configWithPhiSig, philist)
  
odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue, curCov=curCov)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-priorTempered.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig, odemodel)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTempered-",config$kernel,"-",config$seed,".rds"))


muAllDim <- apply(gpode$xsampled, 2:3, mean)
startTheta <- colMeans(gpode$theta)
dotmuAllDim <- apply(gpode$fode, 2:3, mean) 

