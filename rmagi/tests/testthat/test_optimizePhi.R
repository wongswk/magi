#### run with priorTempered phase 1 --------------------------------------------
library(gpds)

config <- list(
  nobs = 33,
  noise = c(0.15,0.15,0.1),
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
  modelName = "Hes1-log",
  startXAtTruth = FALSE,
  startThetaAtTruth = FALSE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  forseTrueMean = FALSE,
  useGPphi1 = FALSE,
  async = TRUE,
  max.epoch = 12,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4],
  phase2 = FALSE
)


config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = log(c(1.438575, 2.037488, 17.90385)),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1logmodelODE(parameters, t(state))))
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
xsim$X3 <- NaN
xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
if(config$async){
  xsim.obs$X1[seq(2,nrow(xsim.obs),by=2)] <- NaN
  xsim.obs$X2[seq(1,nrow(xsim.obs),by=2)] <- NaN
}
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- insertNaN(xsim.obs,config$filllevel)
xsim.obs <- xsim.obs[-nrow(xsim.obs), ]

matplot(xsim$time, xsim[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)
xsim.obs$X2 <- c(xsim.obs$X2[-1], NA)
xsim.obs <- xsim.obs[is.finite(xsim.obs$X1),]

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
xsim.obs$X3 <- (xsim.obs$X2 + xsim.obs$X1)/2
xsim$X3 <- (xsim$X2 + xsim$X1)/2

cursigma <- pram.true$sigma
curphi <- structure(c(2.19417069303522, 63.5419517053017, 0.380334284745747, 
                      33.9780944502718, 0.614883900376933, 68.0653100662423), .Dim = 2:3)

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
xsimInit$X3 <- (xsimInit$X2 + xsimInit$X1)/2
matplot(xsimInit$time, xsimInit[,-1], type="p", pch=2, add=TRUE)
xInit <- data.matrix(xsimInit[,-1])

thetaInit <- c(0.00359214728178393, 0.052904609210373, 0.001, 0.0132093400641409, 
               0.0175094788203553, 0.0205530992047079, 0.001)
sigmaInit <- cursigma


llikXthetaphisigma <- function(xthetaphisigma) {
  xInitial <- matrix(xthetaphisigma[xId], nrow=obsDim[1], ncol=obsDim[2])
  thetaInitial <- xthetaphisigma[thetaId]
  phiInitial <- matrix(xthetaphisigma[phiId], nrow=2)
  sigmaInitial <- xthetaphisigma[sigmaId]
  xthetaphisigmallikRcpp(xInitial, thetaInitial, phiInitial, sigmaInitial,
                         yobs, xsim$time, modelName = "Hes1")
}

phi3optim <- function(xInit, thetaInit, curphi, cursigma){
  fullInit <- c(xInit, thetaInit, curphi, cursigma)
  x3Id <- (length(xInit[, -3]) + 1):length(xInit)
  thetaId <- (max(x3Id)+1):(max(x3Id)+length(thetaInit))
  phi3Id <- (max(thetaId) + length(curphi[,-3]) + 1):(max(thetaId) + length(curphi))
  
  fn <- function(par) {
    fullInit[c(phi3Id)] <- par
    -llikXthetaphisigma( fullInit )$value
  }
  gr <- function(par) {
    fullInit[c(phi3Id)] <- par
    -as.vector(llikXthetaphisigma( fullInit )$grad[c(phi3Id)])
  }
  marlikmap <- optim(curphi[, 3], fn, gr, 
                     method="L-BFGS-B", lower = 0.001, control = list(maxit=1e5))
  matrix(marlikmap$par, nrow=2)
}

xId <- 1:length(data.matrix(xsim[,-1]))
thetaId <- (xId[length(xId)]+1):(xId[length(xId)]+length(pram.true$theta))
phiId <- (thetaId[length(thetaId)]+1):(thetaId[length(thetaId)]+length(pram.true$phi))
sigmaId <- (phiId[length(phiId)]+1):(phiId[length(phiId)]+length(pram.true$sigma))
obsDim <- dim(data.matrix(xsim[,-1]))

phi3list <- phi3optim(xInit, thetaInit, curphi, cursigma)

hes1logmodel <- list(
  fOde=gpds:::hes1logmodelODE,
  fOdeDx=gpds:::hes1logmodelDx,
  fOdeDtheta=gpds:::hes1logmodelDtheta,
  thetaLowerBound=rep(0,7),
  thetaUpperBound=rep(Inf,7)
)


testthat::test_that("optimizePhi can run in R", {
  skip("segfault due to null pointer in hes1logmodel function, passing is wrong")
  gpds:::optimizePhi(yobs, xsim$time, hes1logmodel, cursigma, c(1,1), 
                     xInit, thetaInit, curphi, missingComponentDim=2)
})
