#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = c(0.5, 0.5),
    overrideNoise = TRUE,
    kernel = "finiteDifference2h",
    forceDiagKphi = TRUE,
    forceMean = c("gpsmooth", "truth", "phase8", "zero")[4],
    priorTemperature = c(1, 1),
    seed = 125455454,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 5000,
    burninRatio = 0.20,
    stepSizeFactor = 1,
    filllevel = 2,
    dropoutRate = 0.5,
    modelName = "FN",
    startXAtTruth = TRUE,
    startThetaAtTruth = TRUE,
    startSigmaAtTruth = TRUE,
    sigma_xdot = 0.05
  )
}
config$ndis <- (config$nobs-1)*2^config$filllevel+1

# initialize global parameters, true x, simulated x ----------------------------
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
}
outDir <- paste0(baseDir, "fn-bias/")
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,20,0.05)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
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
curCov <- getCovMphi(config$kernel, xsim, xsim.obs, config=config)
if(config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  for(j in 1:(ncol(xsim)-1)){
    curCov[[j]]$Kinv <- curCov[[j]]$Kinv / config$sigma_xdot^2
    curCov[[j]]$KinvBand <- curCov[[j]]$KinvBand / config$sigma_xdot^2
  }
}
cursigma <- curCov$cursigma
curphi <- curCov$curphi
curCov$cursigma <- NULL
curCov$curphi <- NULL

xtrue.atsim <- sapply(xtrueFunc, function(f) f(xsim$time))
dotxtrue.atsim <- gpds:::fnmodelODE(pram.true$theta, xtrue.atsim)

if(!config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  gpsmoothFuncList <- list()
  gpderivFuncList <- list()
  for(j in 1:(ncol(xsim)-1)){
    ynew <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xtrue$time, 
                         t(curphi[,j]), t(cursigma[j]), kerneltype=config$kernel, deriv = TRUE)
    gpsmoothFuncList[[j]] <- approxfun(xtrue$time, ynew[[1]])
    gpderivFuncList[[j]] <- approxfun(xtrue$time, ynew[[2]])
    plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time),
                  lty = 2, col = j, add = TRUE)
  }
}

if(config$forceMean=="gpsmooth"){
  for(j in 1:(ncol(xsim)-1)){
    curCov[[j]]$mu <- gpsmoothFuncList[[j]](xsim$time)
    curCov[[j]]$dotmu <- gpderivFuncList[[j]](xsim$time)
  }
}else if(config$forceMean=="truth"){
  for(j in 1:(ncol(xsim)-1)){
    curCov[[j]]$mu <- xtrue.atsim[, j]
    curCov[[j]]$dotmu <- dotxtrue.atsim[, j]
  }
  gpsmoothFuncList <- xtrueFunc
  dotxtrue <- gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))
  gpderivFuncList <- lapply(1:ncol(dotxtrue), function(j)
    approxfun(xtrue[, "time"], dotxtrue[, j]))
}else if(config$forceMean=="phase8"){
  phase8 <- readRDS("~/Workspace/DynamicSys/results/phase8-mu.rds")
  curCov[[1]]$mu = phase8[[1]]$mu
  curCov[[1]]$dotmu = phase8[[1]]$dotmu
  curCov[[2]]$mu = phase8[[2]]$mu
  curCov[[2]]$dotmu = phase8[[2]]$dotmu
}else if(config$forceMean=="zero" && config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  gpsmoothFuncList <- list()
  gpderivFuncList <- list()
  for(j in 1:(ncol(xsim)-1)){
    curCov[[j]]$mu = rep(0, nrow(xsim))
    curCov[[j]]$dotmu = rep(0, nrow(xsim))
    gpsmoothFuncList[[j]] <- approxfun(xtrue$time, rep(0, nrow(xtrue)))
    gpderivFuncList[[j]] <- approxfun(xtrue$time, rep(0, nrow(xtrue)))
  }
}

cursigma
curphi
length(curCov)
# MCMC starting value ----------------------------------------------------------
yobs <- data.matrix(xsim[,-1])
xsimInit <- xsim

if(!config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  for(j in 1:(ncol(xsim)-1)){
    nanId <- which(is.na(xsimInit[,j+1]))
    xsimInit[nanId,j+1] <- gpsmoothFuncList[[j]](xsimInit$time[nanId])
  }
  matplot(xsimInit$time, xsimInit[,-1], type="p", pch=2, add=TRUE)
  xInit <- data.matrix(xsimInit[,-1])
  
  thetaInit <- rep(1, length(pram.true$theta))
  
  thetaoptim <- function(xInit, thetaInit, curphi, cursigma){
    curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
      covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                        kerneltype=config$kernel)
      covEach$mu[] <- mean(xsim.obs[,j+1])
      covEach
    })
    fn <- function(par) {
      -xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "FN" )$value
    }
    gr <- function(par) {
      -as.vector(xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "FN" )$grad[-(1:length(xInit))])
    }
    marlikmap <- optim(c(thetaInit), fn, gr, 
                       method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
    thetaInit[] <- marlikmap$par
    list(thetaInit = thetaInit)
  }
  
  thetamle <- thetaoptim(xInit, thetaInit, curphi, cursigma)
  thetaInit <- thetamle$thetaInit
  sigmaInit <- pmax(cursigma, 0.01)
}

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


if(config$forceDiagKphi){
  stopifnot(sum(abs(curCov[[1]]$Kinv)) == sum(abs(diag(curCov[[1]]$Kinv))))
  stopifnot(sum(abs(curCov[[2]]$Kinv)) == sum(abs(diag(curCov[[2]]$Kinv))))
}

# HMC sampler for x, theta, sigma ----------------------------------------------
xId <- 1:length(xInit)
thetaId <- (max(xId)+1):(max(xId)+length(thetaInit))
sigmaId <- (max(thetaId)+1):(max(thetaId)+length(sigmaInit))

if(config$overrideNoise){
  xsim.obs[, -1] <- sapply(xtrueFunc, function(f) f(xsim.obs$time))
  xsim <- insertNaN(xsim.obs,config$filllevel)
  yobs <- data.matrix(xsim[,-1])
}

score_llik <- xthetasigmallikRcpp(
  sapply(xtrueFunc, function(f) f(xsim$time)), theta = pram.true$theta, sigma = pram.true$sigma,
  yobs=yobs, curCov, config$priorTemperature, useBand = TRUE, useMean = TRUE, modelName = config$modelName
)

xthetasigamSingleSampler <- function(xthetasigma, stepSize) {
  observationIndex <- as.numeric(which(apply(yobs, 1, function(yeach) any(is.finite(yeach)))))
  discretizationIndex <- setdiff(1:nrow(xsim), observationIndex)
  sampledIndex <- sample(discretizationIndex, length(discretizationIndex) * (1-config$dropoutRate))
  activeIndex <- sort(c(observationIndex, sampledIndex))
  
  xlasttime <- matrix(xthetasigma[xId], ncol=ncol(yobs))
  activeYobs <- yobs[activeIndex,]
  activeX <- xlasttime[activeIndex,]
  activeStepSize <- matrix(stepSize[xId], ncol=ncol(yobs))[activeIndex,]
  activeCov <- getCovMphi(kernel = config$kernel, xsim = xsim[activeIndex, ], xsim.obs = xsim.obs, config = config)
  for(j in 1:(ncol(xsim)-1)){
    activeCov[[j]]$mu <- curCov[[j]]$mu[activeIndex]
    activeCov[[j]]$dotmu <- curCov[[j]]$dotmu[activeIndex]
  }
  out <- xthetasigmaSample(activeYobs, activeCov, xthetasigma[sigmaId], c(activeX, xthetasigma[thetaId]),
                           c(activeStepSize, stepSize[c(thetaId, sigmaId)]), config$hmcSteps, F, loglikflag = config$loglikflag,
                           priorTemperature = config$priorTemperature, modelName = config$modelName)
  activeXupdated <- matrix(head(out$final, length(activeX)), ncol=ncol(yobs))
  xthistime <- matrix(NA, nrow=nrow(yobs), ncol=ncol(yobs))
  for(j in 1:(ncol(xsim)-1)){
    xthistime[,j] <- approx(xsim$time[activeIndex], activeXupdated[,j], xout=xsim$time)$y
  }
  out$final <- c(xthistime, out$final[-(1:length(activeX))])
  out
}

stepLowInit[sigmaId] <- 0
chainSamplesOut <- chainSampler(config, xthetasigmaInit, xthetasigamSingleSampler, stepLowInit, verbose=TRUE)

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), thetaId],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), xId], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = chainSamplesOut$xth[-(1:burnin), sigmaId, drop=FALSE],
              phi = curphi)
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

configWithPhiSig <- config
if(!config$kernel %in% c("finiteDifference1h", "finiteDifference2h")){
  philist <- lapply(data.frame(round(curphi,3)), function(x) paste(x, collapse = "; "))
  names(philist) <- paste0("phi", 1:length(philist))
  configWithPhiSig <- c(configWithPhiSig, philist)
  configWithPhiSig$score_llik <- score_llik$value
}

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue, curCov=curCov, 
                 gpsmoothFuncList=gpsmoothFuncList, gpderivFuncList=gpderivFuncList)

filename <- paste0(outDir, 
                   config$kernel,"_",
                   "diagKphi-",config$forceDiagKphi,"_",
                   "mu-", config$forceMean,"_",
                   "nobs-", config$nobs,"_",
                   "ndis-", config$ndis,"_",
                   "temperature-xCx-",config$priorTemperature[2], "_",
                   "seed",config$seed)

gpds:::plotPostSamplesFlex(
  paste0(filename,"_gpds.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig, odemodel)

pdf(paste0(filename, "_kernelmat.pdf"))
heatmap_cor(curCov[[1]]$mphi)
title("m-matrix for component V")
heatmap_cor(curCov[[2]]$mphi)
title("m-matrix for component R")
heatmap_cor(curCov[[1]]$Kinv)
title("Kinv-matrix for component V")
heatmap_cor(curCov[[2]]$Kinv)
title("Kinv-matrix for component R")
heatmap_cor(curCov[[1]]$Cinv)
title("Cinv-matrix for component V")
heatmap_cor(curCov[[2]]$Cinv)
title("Cinv-matrix for component R")
dev.off()

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(filename,"_absCI.rds"))
save.image(paste0(filename,"_img.rda"))

