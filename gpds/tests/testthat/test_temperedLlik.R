library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123,
  npostplot = 5,
  loglikflag = "withmean",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 10000,
  burninRatio = 0.1,
  stepSizeFactor = 0.005,
  refitPhiSigma_withmu = FALSE,
  temperatureT = 2000
)

VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=config$noise
)
fn.true <- VRtrue[seq(1,401,by=2),]   #### reference is 201 points
fn.true$time <- seq(0,20,0.1)
fn.sim <- fn.true

set.seed(config$seed)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=config$noise)
tvec.full <- fn.sim$time
fn.sim.all <- fn.sim
fn.sim[-seq(1,nrow(fn.sim), length=config$nobs),] <- NaN
fn.sim.obs <- fn.sim[seq(1,nrow(fn.sim), length=config$nobs),]
tvec.nobs <- fn.sim$time[seq(1,nrow(fn.sim), length=config$nobs)]



foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

if(config$refitPhiSigma_withmu){
  yobsGPsmooth <- (data.matrix(fn.sim[,1:2])-cbind(muV, muR))[!is.nan(fn.sim[,1]),]
}else{
  yobsGPsmooth <- data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2])
}

fn <- function(par) -phisigllikC( par, yobsGPsmooth, 
                                  r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, yobsGPsmooth, 
                                            r.nobs, config$kernel)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
cursigma <- marlikmap$par[5]
curCovV$mu <- as.vector(fn.true[,1])  # pretend these are the means
curCovR$mu <- as.vector(fn.true[,2])

dotmu <- fODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])

cursigma <- mean((data.matrix(fn.sim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(fn.sim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- config$n.iter
xInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)  

config$n.iter <- 99
config$temperatureT <- 1e4
config$stepSizeFactor <- 0.01

stepLowInit <- rep(0.01, 2*nall+3)*config$stepSizeFactor

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag, 
               temperature = config$temperatureT)
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)

llikOut <- xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              chainSamplesOut[[1]][2,])
testthat::expect_equal(llikOut$value/config$temperatureT, chainSamplesOut$lliklist[2])


outsample <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
             chainSamplesOut$xth[1,], chainSamplesOut$stepLow[1, ], 
             config$hmcSteps, T, loglikflag = "band", 
             temperature = 1)

plot.ts(outsample$traj.q[,404])

outsample <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                          xInit, rep(0.001, 2*nall+3), 
                          config$hmcSteps, T, loglikflag = "band", 
                          temperature = 1)
outsample$delta
plot.ts(outsample$traj.H)
plot.ts(outsample$traj.q[,404])

llikOut <- xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              xInit)
llikOut$value
llikOut$grad

#' some breaks early because delta is too large so the proposal must be rejected
#' xthetallik_withmuC, with true mu, likelihood is smaller at biased theta than at true theta
load("/Users/shihaoyang/Workspace/DynamicSys/results/2017-10-29/withmean start mu/first_phase_inference.rda")
xInit1 <- c(muV, muR, startTheta)

#' optimization for whole x does not work, apparently because log posterior is indeed higher at biased point
llikOut <- xthetallikBandApproxC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              xInit)
llikOut$value
llikOut <- xthetallikBandApproxC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              xInit1)
llikOut$value

#' optimization fixing x at posterior mean does not work
llikOut <- xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              xInit)
llikOut$value
llikOut <- xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              xInit1)
llikOut$value

llikOut$grad

fn <- function(par) -xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                                           par)$value
gr <- function(par) -xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                                           par)$grad
marlikmap <- optim(xInit1, fn, gr, method="L-BFGS-B")
marlikmap$par - xInit1
tail(marlikmap$par, 3)
tail(xInit1, 3)

outsample <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                          xInit1, rep(1e-3, 2*nall+3), 
                          config$hmcSteps, T, loglikflag = "withmean", 
                          temperature = 1)
plot.ts(outsample$traj.H)
outsample$acc
outsample$delta
plot.ts(outsample$traj.q[,404])
stepId <- max(which(rowSums(outsample$traj.q) != 0))
plot.ts(outsample$traj.q[stepId,])

llikOut <- xthetallik_withmuC(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma,
                              outsample$traj.q[stepId,])
llikOut$value

plot.ts(outsample$traj.q[,402])

