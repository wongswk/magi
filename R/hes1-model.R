#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
source("R/hes1-helper-functions.R")
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 17,
    noise = c(0.8,0.4,1),
    kernel = "generalMatern",
    seed = 3657260, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 1e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = 1,
    modelName = "Hes1",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = TRUE
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
# config$priorTemperature[2] <- 1e12

# initialize global parameters, true x, simulated x ----------------------------
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                              "-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(0.5, 2, 1),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1modelODE(parameters, t(state))))
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

gpsmoothFuncList <- list()
for(j in 1:3){
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
for(j in 1:3){
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
  with(gpode, gpds:::hes1modelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
names(philist) <- paste0("phi", 1:length(philist))
configWithPhiSig <- c(configWithPhiSig, philist)
  
gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-priorTempered.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTempered-",config$kernel,"-",config$seed,".rds"))


muAllDim <- apply(gpode$xsampled, 2:3, mean)
startTheta <- colMeans(gpode$theta)
dotmuAllDim <- apply(gpode$fode, 2:3, mean) 

# fixing sigma at true value; HMC sampler for x, theta -------------------------
cursigma <- pram.true$sigma
curphi

stepLowInit <- stepLowInit[1:length(c(xInit, thetaInit))]

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, c(xInit, thetaInit), singleSampler, stepLowInit, verbose=TRUE)


gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = matrix(marlikmap$par[-length(marlikmap$par)], 2))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1modelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
configWithPhiSig$sigma <- paste(round(cursigma, 3), collapse = "; ")
philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
names(philist) <- paste0("phi", 1:length(philist))
configWithPhiSig <- c(configWithPhiSig, philist)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-priorTemperedPhase1-trueSigma.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTemperedPhase1-trueSigma-",config$kernel,"-",config$seed,".rds"))


# run with priorTempered phase 2 -------------------------------------------

cursigma <- colMeans((data.matrix(xsim[,-1])-muAllDim)^2, na.rm = TRUE)
cursigma <- sqrt(cursigma)

updatePhi <- curphi

for(j in 1:(ncol(xsim)-1)){
  updatePhiJ <- 
    sapply(as.integer(seq(1, nrow(gpode$xsampled), length=5)), function(postId){
      fn <- function(par) -phillikwithxdotx( par, gpode$xsampled[postId,,j], gpode$fode[postId,,j],
                                             r, signr, config$kernel)
      marlikmap <- optim(rep(10, 2), fn, method="L-BFGS-B", lower = 0.0001,
                         upper = c(Inf, 60*4*2, Inf), hessian = TRUE)
      marlikmap$par  
    })
  updatePhi[,j] <- apply(updatePhiJ, 1, median)
}
updatePhi
curphi
# curphi <- updatePhi

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

for(j in 1:ncol(muAllDim)){
  curCov[[j]]$mu <- muAllDim[,j]
  curCov[[j]]$dotmu <- dotmuAllDim[,j]
}

# refit phi idea doesn't help convergence
# curphi <- matrix(NA, 2, ncol(xsim)-1)
# for(j in 1:(ncol(xsim)-1)){
#   fn <- function(par) -phisigllikC( c(par, cursigma[j]), muAllDim[,j,drop=FALSE],
#                                     r, config$kernel)$value
#   gr <- function(par) -as.vector(phisigllikC( c(par, cursigma[j]), muAllDim[,j,drop=FALSE],
#                                               r, config$kernel)$grad)[1:2]
#   marlikmap <- optim(rep(1, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
#                      upper = c(Inf, 60*4*2))
#   curphi[,j] <- marlikmap$par
# }
# curphi


xInit <- c(muAllDim, startTheta)
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)


gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = matrix(marlikmap$par[-length(marlikmap$par)], 2))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1modelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
configWithPhiSig$sigma <- paste(round(cursigma, 3), collapse = "; ")
philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
names(philist) <- paste0("phi", 1:length(philist))
configWithPhiSig <- c(configWithPhiSig, philist)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-priorTemperedPhase2.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTemperedPhase2-",config$kernel,"-",config$seed,".rds"))

#### last run with true mean --------------------------------------------------

muAllDim <- sapply(xtrueFunc, function(f) f(xsim$time))  # pretend these are the means
dotmuAllDim <- gpds:::hes1modelODE(pram.true$theta, muAllDim) # pretend these are the means for derivatives

for(j in 1:ncol(muAllDim)){
  curCov[[j]]$mu <- muAllDim[,j]
  curCov[[j]]$dotmu <- dotmuAllDim[,j]
}

cursigma <- colMeans((data.matrix(xsim[,-1])-muAllDim)^2, na.rm = TRUE)
cursigma <- sqrt(cursigma)


xInit <- c(muAllDim, pram.true$theta)
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)


gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = matrix(marlikmap$par[-length(marlikmap$par)], 2))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1modelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
configWithPhiSig$sigma <- paste(round(cursigma, 3), collapse = "; ")
philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
names(philist) <- paste0("phi", 1:length(philist))
configWithPhiSig <- c(configWithPhiSig, philist)

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-trueMu.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-trueMu-",config$kernel,"-",config$seed,".rds"))

