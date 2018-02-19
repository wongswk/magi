if(interactive()){
  maxit <- 500
}else{
  maxit <- 10
}

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
  bandsize = 10,
  hmcSteps = 10,
  n.iter = 5000,
  burninRatio = 0.80,
  stepSizeFactor = 1,
  filllevel = 3,
  modelName = "hes1"
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
  fn <- function(par) -phisigloocvllikC( c(par, pram.true$sigma[j]) , data.matrix(xsim.obs[,1+j]), 
                                         r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigloocvllikC( c(par, pram.true$sigma[j]), data.matrix(xsim.obs[,1+j]), 
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


nall <- nrow(xsim)
burnin <- as.integer(config$n.iter*config$burninRatio)
# TODO: add back the pilot idea to solve initial moving to target area problem
# or use tempered likelihood instead
# xInit <- c(vInit, rInit, rep(1,3))
xInit <- c(unlist(lapply(xtrueFunc, function(f) f(xsim$time))), pram.true$theta)
fullInit <- c(xInit, curphi, cursigma)


stepLowInit <- rep(0.000035, length(fullInit))
stepLowInit <- stepLowInit*config$stepSizeFactor

print(llikXthetaphisigma(fullInit))


microbenchmark::microbenchmark(
  singleSamplerLoocvLlik(xInit, rep(1e-5, length(xInit))),
  singleSamplerXthetaphisigma(fullInit, rep(1e-5, length(fullInit)))
)
