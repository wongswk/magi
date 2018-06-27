#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
source("R/hes1-helper-functions.R")
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 33,
    noise = c(0.8,0.4,1.0),
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 2e4,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = 0,
    modelName = "Hes1",
    startXAtTruth = FALSE,
    startThetaAtTruth = FALSE,
    startSigmaAtTruth = TRUE,
    useGPmean = TRUE,
    forseTrueMean = FALSE,
    useGPphi1 = FALSE,
    max.epoch = 12,
    epoch_method = c("mean", "median", "deSolve", "f_x_bar")[1],
    phase2 = FALSE
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
# config$priorTemperature[2] <- config$priorTemperature[1]
# config$priorTemperature <- 1

# initialize global parameters, true x, simulated x ----------------------------
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
  config$seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9 # random seed on cluster
}else{
  baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                              "-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.438575, 2.037488, 17.90385),
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

cursigma <- pram.true$sigma
curphi <- getPhiMarllikFftpriorKnownSigma(cursigma)

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

# fixing sigma at true value; HMC sampler for x, theta -------------------------

cursigma <- pram.true$sigma

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize,
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

if(config$useGPmean){
  for(j in 1:3){
    ydy <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                        t(curphi[,j]), t(cursigma[j]), 
                        kerneltype=config$kernel, deriv = TRUE)
    curCov[[j]]$mu <- as.vector(ydy[[1]])
    curCov[[j]]$dotmu <- as.vector(ydy[[2]])
  }
}

if(config$forseTrueMean){
  dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))
  dotxtrue.atsim <- dotxtrue[xtrue[,"time"] %in% xsim$time,]
  xtrue.atsim <- xtrue[xtrue[,"time"] %in% xsim$time,-1]
  for(j in 1:3){
    curCov[[j]]$mu <- xtrue.atsim[, j]
    curCov[[j]]$dotmu <- dotxtrue.atsim[, j]
  }
}

stepLowInit <- stepLowInit[1:length(c(xInit, thetaInit))]

if(config$useGPphi1){
  curphi[1, 1] <- mean((curCov[[1]]$mu - xsim[,2])^2, na.rm=TRUE)
  curphi[1, 2] <- mean((curCov[[2]]$mu - xsim[,3])^2, na.rm=TRUE)
  curphi[1, 3] <- mean((curCov[[3]]$mu - xsim[,4])^2, na.rm=TRUE)
  
  curCov2 <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize,
                      kerneltype=config$kernel)
    covEach$mu <- curCov[[j]]$mu
    covEach$dotmu <- curCov[[j]]$dotmu
    covEach
  })
  curCov <- curCov2
}


singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, c(xInit, thetaInit), singleSampler, stepLowInit, verbose=TRUE)

## Check sum of square error using numerical solver ------------------------------------
burnin <- as.integer(config$n.iter*config$burninRatio)
mapId <- which.max(chainSamplesOut$lliklist[-(1:burnin)]) + burnin
ttheta <- chainSamplesOut$xth[mapId, (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))]
tx0 <- array(chainSamplesOut$xth[mapId, 1:length(data.matrix(xsim[,-1]))],dim=c(nrow(xsim), ncol(xsim)-1))[1,]
txobs <- array(chainSamplesOut$xth[mapId, 1:length(data.matrix(xsim[,-1]))],dim=c(nrow(xsim), ncol(xsim)-1))[tvec.full %in% tvec.nobs,]

xtrue2 <- deSolve::ode(y = tx0, times = times, func = modelODE, parms = ttheta)

matplot(xtrue[, "time"], xtrue2[, -1], type="l", lty=1)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=3, add=T)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
xtrue.obs <- xtrue[xtrue[,"time"] %in% xsim.obs$time,-1]
xsampler.obs <- xtrue2[xtrue2[,"time"] %in% xsim.obs$time,-1]

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue, curCov=curCov)

rmseTrue <- sqrt(apply((xtrue.obs[,1:2] - xsim.obs[,2:3])^2,2,mean))
rmseWholeGpode <- sqrt(apply((txobs[,1:2] - xsim.obs[,2:3])^2,2, mean))
rmseOdeGpode <- sqrt(apply((xsampler.obs[,1:2] - xsim.obs[,2:3])^2,2,mean))

rmselist <- list(
  rmseTrue = paste0(round(rmseTrue, 3), collapse = "; "),
  rmseWholeGpode = paste0(round(rmseWholeGpode, 3), collapse = "; "),
  rmseOdeGpode = paste0(round(rmseOdeGpode, 3), collapse = "; ")
)

#### end SSE check

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = curphi)
sampleId <- dim(gpode$xsampled)[1]
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
  paste0(outDir, config$kernel,"-",config$seed,"-",date(),"-hes1-fully-observed-trueSigma.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig, odemodel)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-hes1-fully-observed-trueSigma-",config$kernel,"-",config$seed,".rds"))
save.image(paste0(
  outDir,
  config$loglikflag,"-hes1-fully-observed-trueSigma-",config$kernel,"-",config$seed,".rda"))

pdf(paste0(outDir, config$kernel,"-",config$seed,"-",date(),"-hes1-fully-observed-trueSigma-addon.pdf"), 
    width = 8, height = 8)
layout(matrix(1:4, 2))
matplot(xtrue$time, xtrue[,-1], type="l", lty=1, main="raw data")
matplot(xsim$time, xsim[,-1], type="p", pch=20, add=TRUE)

npostplot=50
id.plot <- seq(1,nrow(gpode$theta),length=npostplot)
id.plot <- unique(as.integer(id.plot))

for(j in 1:(ncol(xsim)-1)){
  matplot(xtrue$time, cbind(xtrue[,j+1], dotxtrue[,j]), type="l", lty=1, col=c(2,1),
          ylab=paste0("component-",j), main="full posterior")
  mtext(paste("component", c("P", "M", "H")[j]))
  points(xsim$time, xsim[,j+1], col=2)
  matplot(xsim$time, t(gpode$xsampled[id.plot,,j]), col="skyblue",add=TRUE, type="b",lty=1, pch=20)
  matplot(xsim$time, t(gpode$fode[id.plot,,j]), col="grey",add=TRUE, type="b",lty=1, pch=20)
  
  if(!is.null(odemodel) && !is.null(odemodel$curCov)){
    lines(xsim$time, odemodel$curCov[[j]]$mu, col="forestgreen", lwd=2)
    lines(xsim$time, odemodel$curCov[[j]]$dotmu, col="darkgreen", lwd=2)
  }
}

layout(matrix(1:(2*length(pram.true$theta)),ncol=2,byrow = TRUE))
par(mar=c(2,1,1.5,1))
for(i in 1:length(pram.true$theta)){
  hist(gpode$theta[,i], main=letters[i], xlab=NA, ylab=NA)
  abline(v=pram.true$theta[i], col=2)
  abline(v=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
  plot.ts(gpode$theta[,i], main=letters[i], xlab=NA, ylab=NA)
  abline(h=pram.true$theta[i], col=2)
  abline(h=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
}

dev.off()

# phase 2 ----------------------------------------------------------------------
if(config$phase2){
muAllDim <- apply(gpode$xsampled, 2:3, mean)
dotmuAllDim <- apply(gpode$fode, 2:3, mean) 

for(j in 1:ncol(muAllDim)){
  curCov[[j]]$mu <- muAllDim[,j]
  curCov[[j]]$dotmu <- dotmuAllDim[,j]
}

xInit <- muAllDim
thetaInit <- colMeans(gpode$theta)
stepLowInit <- chainSamplesOut$stepLow
config$burninRatio <- 0.5
config$phase2 <- TRUE

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, c(xInit, thetaInit), singleSampler, stepLowInit, verbose=TRUE)

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue, curCov=curCov)

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = curphi)
sampleId <- dim(gpode$xsampled)[1]
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
  paste0(outDir, config$kernel,"-",config$seed,"-",date(),"-hes1-fully-observed-trueSigma-phase2.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig, odemodel)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-hes1-fully-observed-trueSigma-phase2-",config$kernel,"-",config$seed,".rds"))
save.image(paste0(
  outDir,
  config$loglikflag,"-hes1-fully-observed-trueSigma-phase2-",config$kernel,"-",config$seed,".rda"))

pdf(paste0(outDir, config$kernel,"-",config$seed,"-",date(),"-hes1-fully-observed-trueSigma-phase2-addon.pdf"), 
    width = 8, height = 8)
layout(matrix(1:4, 2))
matplot(xtrue$time, xtrue[,-1], type="l", lty=1, main="raw data")
matplot(xsim$time, xsim[,-1], type="p", pch=20, add=TRUE)

npostplot=50
id.plot <- seq(1,nrow(gpode$theta),length=npostplot)
id.plot <- unique(as.integer(id.plot))

for(j in 1:(ncol(xsim)-1)){
  matplot(xtrue$time, cbind(xtrue[,j+1], dotxtrue[,j]), type="l", lty=1, col=c(2,1),
          ylab=paste0("component-",j), main="full posterior")
  mtext(paste("component", c("P", "M", "H")[j]))
  points(xsim$time, xsim[,j+1], col=2)
  matplot(xsim$time, t(gpode$xsampled[id.plot,,j]), col="skyblue",add=TRUE, type="b",lty=1, pch=20)
  matplot(xsim$time, t(gpode$fode[id.plot,,j]), col="grey",add=TRUE, type="b",lty=1, pch=20)
  
  if(!is.null(odemodel) && !is.null(odemodel$curCov)){
    lines(xsim$time, odemodel$curCov[[j]]$mu, col="forestgreen", lwd=2)
    lines(xsim$time, odemodel$curCov[[j]]$dotmu, col="darkgreen", lwd=2)
  }
}

layout(matrix(1:(2*length(pram.true$theta)),ncol=2,byrow = TRUE))
par(mar=c(2,1,1.5,1))
for(i in 1:length(pram.true$theta)){
  hist(gpode$theta[,i], main=letters[i], xlab=NA, ylab=NA)
  abline(v=pram.true$theta[i], col=2)
  abline(v=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
  plot.ts(gpode$theta[,i], main=letters[i], xlab=NA, ylab=NA)
  abline(h=pram.true$theta[i], col=2)
  abline(h=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
}

dev.off()
}

# epoch ----------------------------------------------------------------------------------
config$burninRatio <- 0.2
config$n.iter <- 5000
epoch_sequence <- 3:config$max.epoch
epoch_sequence <- epoch_sequence[epoch_sequence > 3]
thismodelODE <- gpds:::hes1modelODE

for(epoch in epoch_sequence){
  if(config$epoch_method == "mean"){
    muAllDim <- apply(gpode$xsampled, 2:3, mean)
    dotmuAllDim <- apply(gpode$fode, 2:3, mean) 
  }else if(config$epoch_method == "f_x_bar"){
    muAllDim <- apply(gpode$xsampled, 2:3, mean)
    thetaInit <- apply(gpode$theta, 2, mean)
    dotmuAllDim <- thismodelODE(thetaInit, muAllDim)
  }else if(config$epoch_method == "median"){
    muAllDim <- apply(gpode$xsampled, 2:3, median)
    dotmuAllDim <- apply(gpode$fode, 2:3, median) 
  }else if(config$epoch_method == "deSolve"){
    ttheta <- colMeans(gpode$theta)
    tx0 <- colMeans(gpode$xsampled[,1,])
    xdesolvePM <- deSolve::ode(y = tx0, times = xsim$time, func = modelODE, parms = ttheta)
    muAllDim <- as.matrix(xdesolvePM[,-1])
    dotmuAllDim <- thismodelODE(ttheta, muAllDim)
  }
  
  for(j in 1:ncol(muAllDim)){
    curCov[[j]]$mu <- muAllDim[,j]
    curCov[[j]]$dotmu <- dotmuAllDim[,j]
  }
  
  xInit <- muAllDim
  thetaInit <- apply(gpode$theta, 2, mean)
  stepLowInit <- chainSamplesOut$stepLow
  
  singleSampler <- function(xthetaValues, stepSize) 
    xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
                 xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
                 priorTemperature = config$priorTemperature, modelName = "Hes1")
  chainSamplesOut <- chainSampler(config, c(xInit, thetaInit), singleSampler, stepLowInit, verbose=TRUE)
  
  odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue, curCov=curCov)
  
  burnin <- as.integer(config$n.iter*config$burninRatio)
  gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
                xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                               dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
                lglik=chainSamplesOut$lliklist[-(1:burnin)],
                sigma = cursigma,
                phi = curphi)
  sampleId <- dim(gpode$xsampled)[1]
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
    paste0(outDir, config$kernel,"-",config$seed,"-",date(),"-hes1-fully-observed-trueSigma-phase",epoch,".pdf"), 
    xtrue, dotxtrue, xsim, gpode, pram.true, configWithPhiSig, odemodel)
  
  absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
  absCI <- rbind(absCI, mean=colMeans(gpode$theta))
  absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))
  
  saveRDS(absCI, paste0(
    outDir,
    config$loglikflag,"-hes1-fully-observed-trueSigma-phase",epoch,"-",config$kernel,"-",config$seed,".rds"))
  save.image(paste0(
    outDir,
    config$loglikflag,"-hes1-fully-observed-trueSigma-phase",epoch,"-",config$kernel,"-",config$seed,".rda"))
  
}
