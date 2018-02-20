#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
library(parallel)

config <- list(
  nobs = 11,
  noise = c(4,1,8)*0.2,
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2000,
  burninRatio = 0.80,
  stepSizeFactor = 1,
  filllevel = 3,
  modelName = "Hes1",
  pilotSize = 40,
  pilotRatio = 1.0
)


config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/"  
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

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                    r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                              r.nobs, config$kernel)$grad)
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}
cursigmaMarginalLikelihood <- cursigma
curphiMarginalLikelihood <- curphi

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), 
                                         r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), 
                                                   r.nobs, config$kernel)$grad)
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}
cursigmaLoocvLlik <- cursigma
curphiLoocvLlik <- curphi

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -phisigloocvmseC( c(par, pram.true$sigma[j]) , data.matrix(xsim.obs[,1+j]), 
                                         r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigloocvmseC( c(par, pram.true$sigma[j]), data.matrix(xsim.obs[,1+j]), 
                                                   r.nobs, config$kernel)$grad)[1:2]
  marlikmap <- optim(rep(100, 2), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- pram.true$sigma[j]
  curphi[,j] <- marlikmap$par[1:2]
}
cursigmaLoocvMse <- cursigma
curphiLoocvMse <- curphi

rm(cursigma, curphi)


curCovLoocvMse <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphiLoocvMse[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})
curCovLoocvLlik <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphiLoocvLlik[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})
curCovMarginalLikelihood <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphiMarginalLikelihood[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})


singleSamplerLoocvLlik <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCovLoocvLlik, cursigmaLoocvLlik, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
singleSamplerLoocvMse <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCovLoocvMse, cursigmaLoocvMse, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
singleSamplerMarginalLikelihood <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCovMarginalLikelihood, cursigmaMarginalLikelihood, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")

yobs <- data.matrix(xsim[,-1])
obsDim <- dim(yobs)
xId <- 1:length(yobs)
thetaId <- (xId[length(xId)]+1):(xId[length(xId)]+length(pram.true$theta))
phiId <- (thetaId[length(thetaId)]+1):(thetaId[length(thetaId)]+length(pram.true$phi))
sigmaId <- (phiId[length(phiId)]+1):(phiId[length(phiId)]+length(pram.true$sigma))

singleSamplerXthetaphisigma <- function(xthetaphisigma, stepSize) {
  print(Sys.time())
  xInitial <- matrix(xthetaphisigma[xId], nrow=obsDim[1], ncol=obsDim[2])
  thetaInitial <- xthetaphisigma[thetaId]
  phiInitial <- matrix(xthetaphisigma[phiId], nrow=2)
  sigmaInitial <- xthetaphisigma[sigmaId]
  xthetaphisigmaSample(xInitial, thetaInitial, phiInitial, sigmaInitial,
                       yobs, xsim$time, stepSize, modelName = "Hes1", nsteps = config$hmcSteps)
}

llikXthetaphisigma <- function(xthetaphisigma) {
  xInitial <- matrix(xthetaphisigma[xId], nrow=obsDim[1], ncol=obsDim[2])
  thetaInitial <- xthetaphisigma[thetaId]
  phiInitial <- matrix(xthetaphisigma[phiId], nrow=2)
  sigmaInitial <- xthetaphisigma[sigmaId]
  xthetaphisigmallikRcpp(xInitial, thetaInitial, phiInitial, sigmaInitial,
                         yobs, xsim$time, modelName = "Hes1")
}


rbind(curphiLoocvMse, cursigmaLoocvMse, 
      curphiLoocvLlik, cursigmaLoocvLlik,
      curphiMarginalLikelihood, cursigmaMarginalLikelihood)





singleSampler <- singleSamplerXthetaphisigma
cursigma <- apply(cbind(cursigmaLoocvLlik, cursigmaLoocvMse, cursigmaMarginalLikelihood), 1, median)
curphi <- apply(abind::abind(curphiLoocvLlik, curphiLoocvMse, curphiMarginalLikelihood, along=3), 1:2, median)

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    penalty <- dnorm(par[2], max(xsim.obs$time)/2, max(xsim.obs$time)/6, log=TRUE)
    -(marlik$value + loocvlik$value) # + penalty
  }
  gr <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad + loocvlik$grad)
    # grad[2] <- grad[2] + (par[2] - max(xsim.obs$time)/2) / (max(xsim.obs$time)/6)^2
    grad
  }
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}

cursigma
curphi

# gpsmooth function ------------------------------------------------------------
gpsmoothFuncList <- list()
for(j in 1:3){
  ynew <- getMeanCurve(xsim.obs$time, xsim.obs[,j+1], xsim$time, 
                       t(curphi[,j]), t(cursigma[j]), kerneltype=config$kernel)
  gpsmoothFuncList[[j]] <- approxfun(xsim$time, ynew)
  plot.function(gpsmoothFuncList[[j]], from = min(xsim$time), to = max(xsim$time),
                lty = 2, col = j, add = TRUE)
}

# set thetaInit by optim -------------------------------------------------------
thetaoptim <- function(xInit, thetaInit, curphi, cursigma){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    -xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "Hes1" )$value
  }
  gr <- function(par) {
    -as.vector(xthetallikRcpp( yobs, curCov, cursigma, c(xInit, par), "Hes1" )$grad[-(1:length(xInit))])
  }
  marlikmap <- optim(c(thetaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  thetaInit[] <- marlikmap$par
  list(thetaInit = thetaInit)
}

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

# visualization utilities ------------------------------------------------------
plotXinit <- function(...){
  xInitList <- list(...)
  matplot(xsim.obs$time, xsim.obs[,-1], col=1:3, pch=20)
  for(j in 1:length(xInitList)){
    matplot(xsim$time, xInitList[[j]], col=1:3, add = TRUE, 
            type="p", lwd=1, lty=2, pch=j)
  }
}


# optimize phi -----------------------------------------------------------------

llikXthetaphisigma(c(xInit, thetaInit, curphi, cursigma))

# ususal R optim
fulloptim <- function(xInit, thetaInit, curphi, cursigma){
  fn <- function(par) -llikXthetaphisigma( par )$value
  gr <- function(par) -as.vector(llikXthetaphisigma( par )$grad)
  marlikmap <- optim(c(xInit, thetaInit, curphi, cursigma), 
                     fn, gr, method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  
  xInit[] <- marlikmap$par[xId]
  thetaInit[] <- marlikmap$par[thetaId]
  curphi[] <- marlikmap$par[phiId]
  cursigma[] <- marlikmap$par[sigmaId]
  list(xInit = xInit, 
       thetaInit = thetaInit, 
       curphi = curphi, 
       cursigma = cursigma)
}
# fullmle <- fulloptim(data.matrix(xsimInit[,-1]), thetaInit, curphi, cursigma)  # unstable, could take long

# matplot(xsim$time, fullmle$xInit, lty=4, col=1:3, add = TRUE, type="p", lwd=3, pch=4)

# R optim with iterative optimization
phisigmaoptim <- function(xInit, thetaInit, curphi, cursigma){
  fullInit <- c(xInit, thetaInit, curphi, cursigma)
  fn <- function(par) {
    fullInit[c(phiId, sigmaId)] <- par
    -llikXthetaphisigma( fullInit )$value
  }
  gr <- function(par) {
    fullInit[c(phiId, sigmaId)] <- par
    -as.vector(llikXthetaphisigma( fullInit )$grad[c(phiId, sigmaId)])
  }
  marlikmap <- optim(c(curphi, cursigma), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  curphi[] <- marlikmap$par[1:length(curphi)]
  cursigma[] <- marlikmap$par[-(1:length(curphi))]
  list(curphi = curphi,
       cursigma = cursigma)
}

xthetaoptim <- function(xInit, thetaInit, curphi, cursigma, priorTemperature = rep(1,2)){
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  fn <- function(par) {
    -xthetallikRcpp( yobs, curCov, cursigma, par, "Hes1", priorTemperatureInput=priorTemperature )$value
  }
  gr <- function(par) {
    -as.vector(xthetallikRcpp( yobs, curCov, cursigma, par, "Hes1", priorTemperatureInput=priorTemperature )$grad)
  }
  marlikmap <- optim(c(xInit, thetaInit), fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  xInit[] <- marlikmap$par[1:length(xInit)]
  thetaInit[] <- marlikmap$par[-(1:length(xInit))]
  list(xInit = xInit,
       thetaInit = thetaInit)
}

fullInit <- list(xInit=xInit, 
                 thetaInit=thetaInit, 
                 curphi=curphi, 
                 cursigma=cursigma)

fullInit[-1]
plotXinit(fullInit$xInit)

heatUpBeta <- matrix(c(1, 1), nrow = 2, ncol = 3)
priorTemperature <- c(1, 1e2)
xthetamle <- with(fullInit, xthetaoptim(xInit, thetaInit, heatUpBeta*curphi, cursigma, priorTemperature))
plotXinit(fullInit$xInit, xthetamle$xInit)

fullInit[names(xthetamle)] <- xthetamle
priorTemperature <- c(1, 1)
xthetamle <- with(fullInit, xthetaoptim(xInit, thetaInit, heatUpBeta*curphi, cursigma, priorTemperature))
plotXinit(fullInit$xInit, xthetamle$xInit)

fullInit[names(xthetamle)] <- xthetamle
thetamle <- with(fullInit, thetaoptim(xInit, thetaInit, curphi, cursigma))
thetamle$thetaInit - fullInit$thetaInit

phisigmamle <- with(fullInit, phisigmaoptim(xInit, thetaInit, curphi, cursigma))
phisigmamle
fullInit[names(phisigmamle)] <- phisigmamle


phisigmamle
matplot(xsim$time, xthetamle$xInit, col=1:3, add = TRUE, type="p", lwd=3, pch=5)

llikXthetaphisigma(unlist(fullInit))



# stochastic gradient descend for optimization task ----------------------------
# initllikOld <- list()
# fullInit <- marlikmap$par
# 
# learningRate <- rep(1e-4, length(fullInit))
# learningRate[xId] <- 0
# learningRate[thetaId] <- 0
# learningRate[phiId] <- 1e-4
# learningRate[sigmaId] <- 1e-5
# 
# for(sgdIt in 1:1e6){
#   initllik <- llikXthetaphisigma(fullInit)
#   print(initllik$value - initllikOld$value)
#   (initllik$value - initllikOld$value)/abs(initllik$value)
#   initllikOld <- initllik
#   summary(initllik$grad[learningRate>0])
#   fullInit <- fullInit + initllik$grad*learningRate
# }
# 
# xthetaphisigma <- fullInit
# xInitial <- matrix(xthetaphisigma[xId], nrow=obsDim[1], ncol=obsDim[2])
# thetaInitial <- xthetaphisigma[thetaId]
# phiInitial <- matrix(xthetaphisigma[phiId], nrow=2)
# sigmaInitial <- xthetaphisigma[sigmaId]
# 
# matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
# matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
# matplot(xsim$time, xInitial, type="p", col=1:(ncol(xsim)-1), pch=3, add = TRUE)



# E-M like iterative optim -----------------------------------------------------
# 1. conditional on phi-sigma, sample x-theta
# 2. approximate E-step
# 3. optim M-step


# x theta posterior sampler ----------------------------------------------------

xthetaInit <- unlist(fullInit[c("xInit", "thetaInit")])
stepLowInit <- rep(config$stepSize, length(xthetaInit))

curphi <- fullInit$curphi
cursigma <- fullInit$cursigma

curCov <- lapply(1:(ncol(xsim)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, xthetaInit, singleSampler, stepLowInit, verbose=TRUE)


nall <- nrow(xsim)
burnin <- as.integer(config$n.iter*config$burninRatio)


gpode <- list(theta=chainSamplesOut$xth[-(1:burnin), (length(data.matrix(xsim[,-1]))+1):(ncol(chainSamplesOut$xth))],
              xsampled=array(chainSamplesOut$xth[-(1:burnin), 1:length(data.matrix(xsim[,-1]))], 
                             dim=c(config$n.iter-burnin, nall, ncol(xsim)-1)),
              lglik=chainSamplesOut$lliklist[-(1:burnin)],
              sigma = cursigma,
              phi = matrix(marlikmap$par[-length(marlikmap$par)], 2))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1modelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1modelODE(pram.true$theta, data.matrix(xtrue[,-1]))

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$kernel,"-",config$seed,"-priorTempered.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config)

absCI <- apply(gpode$theta, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$theta))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$theta &  pram.true$theta < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTempered-",config$kernel,"-",config$seed,".rds"))


muAllDim <- apply(gpode$xsampled, 2:3, mean)
startTheta <- colMeans(gpode$theta)
dotmuAllDim <- apply(gpode$fode, 2:3, mean) 
