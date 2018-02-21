#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
source("R/hes1-helper-functions.R")
config <- list(
  nobs = 11,
  noise = c(4,1,8)*0.2,
  kernel = "generalMatern",
  seed = 3657260, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 1e4,
  burninRatio = 0.50,
  stepSizeFactor = 1,
  filllevel = 3,
  modelName = "Hes1"
)


config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
# config$priorTemperature[2] <- 1e12

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

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
  
  fn <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                     max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
    -(marlik$value + loocvlik$value + penalty)
  }
  gr <- function(par) {
    marlik <- phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    loocvlik <- phisigloocvllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad + loocvlik$grad)
    grad[2] <- grad[2] + (par[2] - priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
    grad
  }
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}

cursigma
curphi

eval(phiAll3methodsExpr)


phi2candidates <- exp(seq(0, 6, length=100))
phi2candidatesCheckedMarllik <- sapply(phi2candidates, checkPhi2FitMarllik, j=2)
phi2candidatesCheckedMarllikFuncList <-
  apply(phi2candidatesCheckedMarllik, 2, function(phisig){
    getGPsmoothFunc(phisig[2:3], phisig[4], j=3, showplot = FALSE)  
  })

phi2candidatesCheckedLoocvllik <- sapply(phi2candidates, checkPhi2FitLoocvllik, j=2)
phi2candidatesCheckedLoocvllikFuncList <-
  apply(phi2candidatesCheckedLoocvllik, 2, function(phisig){
    getGPsmoothFunc(phisig[2:3], phisig[4], j=3, showplot = FALSE)  
  })

phi2candidatesCheckedLoocvmse <- sapply(phi2candidates, checkPhi2FitLoocvmse, j=2)
phi2candidatesCheckedLoocvmseFuncList <-
  apply(phi2candidatesCheckedLoocvmse, 2, function(phisig){
    getGPsmoothFunc(phisig[2:3], phisig[4], j=3, showplot = FALSE)  
  })

pdf("fitting based on only observations.pdf", height = 16, width = 12)
layout(matrix(1:12, nrow=4))
plot(phi2candidates, phi2candidatesCheckedMarllik[1,], main="Marllik likelihood value")
plot(xsim.obs$time, xsim.obs[,4], pch=20, col=3, main="Marllik fitted curve")
sapply(1:length(phi2candidatesCheckedMarllikFuncList), function(id){
  f <- phi2candidatesCheckedMarllikFuncList[[id]]
  plot.function(f, from = min(xsim$time), to = max(xsim$time),
                lty = 1, col = paste0("gray",id), add = TRUE)
})
points(xsim.obs$time, xsim.obs[,4], pch=20, col=3)
plot(phi2candidates, phi2candidatesCheckedMarllik[2,], main="Marllik phi1")
plot(phi2candidates, phi2candidatesCheckedMarllik[4,], main="Marllik sigma")

plot(phi2candidates, phi2candidatesCheckedLoocvllik[1,], main="Loocvllik likelihood value")
plot(xsim.obs$time, xsim.obs[,4], pch=20, col=3, main="Loocvllik fitted curve")
sapply(1:length(phi2candidatesCheckedLoocvllikFuncList), function(id){
  f <- phi2candidatesCheckedLoocvllikFuncList[[id]]
  plot.function(f, from = min(xsim$time), to = max(xsim$time),
                lty = 1, col = paste0("gray",id), add = TRUE)
})
points(xsim.obs$time, xsim.obs[,4], pch=20, col=3)
plot(phi2candidates, phi2candidatesCheckedLoocvllik[2,], main="Loocvllik phi1")
plot(phi2candidates, phi2candidatesCheckedLoocvllik[4,], main="Loocvllik sigma")

plot(phi2candidates, phi2candidatesCheckedLoocvmse[1,], main="Loocvmse likelihood value")
plot(xsim.obs$time, xsim.obs[,4], pch=20, col=3, main="Loocvmse fitted curve")
sapply(1:length(phi2candidatesCheckedLoocvmseFuncList), function(id){
  f <- phi2candidatesCheckedLoocvmseFuncList[[id]]
  plot.function(f, from = min(xsim$time), to = max(xsim$time),
                lty = 1, col = paste0("gray",id), add = TRUE)
})
points(xsim.obs$time, xsim.obs[,4], pch=20, col=3)
plot(phi2candidates, phi2candidatesCheckedLoocvmse[2,], main="Loocvmse phi1")
plot(phi2candidates, phi2candidatesCheckedLoocvmse[4,], main="Loocvmse sigma")
dev.off()


curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

nall <- nrow(xsim)
burnin <- as.integer(config$n.iter*config$burninRatio)

xInit <- c(unlist(lapply(xtrueFunc, function(f) f(xsim$time))), pram.true$theta)
stepLowInit <- rep(0.000035, (ncol(xsim)-1)*nall+length(pram.true$theta))
stepLowInit <- stepLowInit*config$stepSizeFactor

# TODO: add back the pilot idea to solve initial moving to target area problem
# or use tempered likelihood instead
# xInit <- c(vInit, rInit, rep(1,3))

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,-1]), curCov, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature, modelName = "Hes1")
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)


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

configWithPhiSig <- config
configWithPhiSig$sigma <- paste(round(cursigma, 3), collapse = "; ")
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
