library(magi)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123,
  npostplot = 5,
  loglikflag = "band",
  bandsize = 20,
  hmcSteps = 20,
  n.iter = 2e2,
  burninRatio = 0.1,
  stepSizeFactor = 1
)

VRtrue <- read.csv(system.file("testdata/FN.csv", package="magi"))
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

priorFactor <- magi:::calcFrequencyBasedPrior(fn.sim.obs[,1])
priorFactor2 <- magi:::getFrequencyBasedPrior(fn.sim.obs[,1])

testthat::test_that("c++ calcFrequencyBasedPrior correct", {
  testthat::expect_true(all(abs(priorFactor - priorFactor2) < 1e-3))
})

testthat::test_that("c++ gpsmooth correct", {
  r.nobs <- abs(outer(tvec.nobs, t(tvec.nobs),'-')[,1,])
  yobs1 <- data.matrix(fn.sim.obs[,1,drop=FALSE])
  outputc <- magi:::gpsmooth(yobs1,
                             r.nobs,
                             config$kernel)
  
  fn <- function(par) {
    marlik <- magi:::phisigllikC( par, yobs1, r.nobs, config$kernel)
    -marlik$value
  }
  gr <- function(par) {
    marlik <- magi:::phisigllikC( par, yobs1, r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad)
    grad
  }

  fn(outputc)
  testthat::expect_true(all(abs(gr(outputc)) < 1e-3))
})

testthat::test_that("c++ gpsmooth correct dim2", {
  r.nobs <- abs(outer(tvec.nobs, t(tvec.nobs),'-')[,1,])
  yobs1 <- data.matrix(fn.sim.obs[,1:2,drop=FALSE])
  outputc <- magi:::gpsmooth(yobs1,
                             r.nobs,
                             config$kernel)
  
  fn <- function(par) {
    marlik <- magi:::phisigllikC( par, yobs1, r.nobs, config$kernel)
    -marlik$value
  }
  gr <- function(par) {
    marlik <- magi:::phisigllikC( par, yobs1, r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad)
    grad
  }
  
  fn(outputc)
  testthat::expect_true(all(abs(gr(outputc)) < 1e-3))
})

testthat::test_that("c++ gpsmooth correct with fft prior", {
  r.nobs <- abs(outer(tvec.nobs, t(tvec.nobs),'-')[,1,])
  yobs1 <- data.matrix(fn.sim.obs[,1,drop=FALSE])
  outputc <- magi:::gpsmooth(yobs1,
                             r.nobs,
                             config$kernel,
                             -1,
                             TRUE)
  xsim.obs <- fn.sim.obs[,c("time", "Vtrue", "Rtrue")]
  j=1
  fn <- function(par) {
    marlik <- magi:::phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor[1], 
                     max(xsim.obs$time)*priorFactor[2], log=TRUE)
    # penalty <- dgamma(par[2], alphaRate, betaRate/max(xsim.obs$time), log=TRUE)
    # penalty <- 0
    -(marlik$value + penalty)
  }
  gr <- function(par) {
    marlik <- magi:::phisigllikC( par, data.matrix(xsim.obs[,1+j]), r.nobs, config$kernel)
    grad <- -as.vector(marlik$grad)
    penalty <- (par[2] - max(xsim.obs$time)*priorFactor[1]) / (max(xsim.obs$time)*priorFactor[2])^2
    # penalty <- ((alphaRate-1)/par[2] - betaRate/max(xsim.obs$time))
    # penalty <- 0
    grad[2] <- grad[2] + penalty
    grad
  }
  
  fn(outputc)
  testthat::expect_true(all(abs(gr(outputc)) < 1e-4))
})

testthat::test_that("c++ gpsmooth correct with fixed sigma fft prior", {
  r.nobs <- abs(outer(tvec.nobs, t(tvec.nobs),'-')[,1,])
  yobs1 <- data.matrix(fn.sim.obs[,1,drop=FALSE])
  outputc <- magi:::gpsmooth(yobs1,
                             r.nobs,
                             config$kernel,
                             0.1,
                             FALSE)
  testthat::expect_equal(length(outputc), 2)
  
  outputc <- magi:::gpsmooth(yobs1,
                             r.nobs,
                             config$kernel,
                             0.1,
                             TRUE)
  testthat::expect_equal(length(outputc), 2)
  
})

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

marlikmap <- list(par=c(2.314334, 1.346233, 0.622316, 2.451729, 0.084745))
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
cursigma <- marlikmap$par[5]
curCovV$mu <- as.vector(fn.true[,1])  # pretend these are the means
curCovR$mu <- as.vector(fn.true[,2])
curCovV$tvecCovInput <- tvec.full
curCovR$tvecCovInput <- tvec.full

dotmu <- magi:::fODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])

nall <- nrow(fn.sim)
burnin <- as.integer(config$n.iter*config$burninRatio)

xInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

singleSampler <- function(xthetaValues, stepSize)
  magi:::xthetaSample(data.matrix(fn.sim[,1:2]), list(curCovV, curCovR), cursigma,
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag)

testthat::test_that("chainSampler can run without error",{
  chainSamplesOut <- magi:::chainSampler(config, xInit, singleSampler, stepLowInit, verbose=FALSE)
  testthat::expect_equal(length(chainSamplesOut$lliklist), config$n.iter)
})

fnmodel <- list(
  fOde=magi:::fnmodelODE,
  fOdeDx=magi:::fnmodelDx,
  fOdeDtheta=magi:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf)
)

testthat::test_that("chainSamplerRcpp can run without error",{
  xthetasigmaInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc, c(cursigma, cursigma))
  stepLowXthetasigmaInit <- c(rep(0.00035, 2*nall+3)*config$stepSizeFactor, 0, 0)

  out <- magi:::chainSamplerRcpp(
    yobs = data.matrix(fn.sim[,1:2]),
    covAllDimInput = list(curCovV, curCovR),
    nstepsInput = config$hmcSteps,
    loglikflagInput = config$loglikflag,
    priorTemperatureInput = c(1, 1),
    sigmaSizeInput = 2,
    modelInput = fnmodel,
    niterInput = config$n.iter,
    burninRatioInput = config$burninRatio,
    xthetasigmaInit = xthetasigmaInit,
    stepLowInit = stepLowXthetasigmaInit,
    verbose = TRUE
  )
  testthat::expect_equal(length(out$lliklist), config$n.iter)
  
  xthetasigmaInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc, c(cursigma))
  stepLowXthetasigmaInit <- c(rep(0.00035, 2*nall+3)*config$stepSizeFactor, 0)
  out <- magi:::chainSamplerRcpp(
    yobs = data.matrix(fn.sim[,1:2]),
    covAllDimInput = list(curCovV, curCovR),
    nstepsInput = config$hmcSteps,
    loglikflagInput = config$loglikflag,
    priorTemperatureInput = c(1, 1),
    sigmaSizeInput = 1,
    modelInput = fnmodel,
    niterInput = config$n.iter,
    burninRatioInput = config$burninRatio,
    xthetasigmaInit = xthetasigmaInit,
    stepLowInit = stepLowXthetasigmaInit,
    verbose = TRUE
  )
  testthat::expect_equal(length(out$lliklist), config$n.iter)

})

testthat::test_that("optimizeThetaInit can run without error",{
   out <- magi:::optimizeThetaInitRcpp(
     yobs = data.matrix(fn.sim[,1:2]),
     odeModel = fnmodel,
     covAllDimInput = list(curCovV, curCovR),
     sigmaAllDimensionsInput = c(cursigma, cursigma),
     priorTemperatureInput = c(1,1),
     xInitInput = cbind(fn.true$Vtrue, fn.true$Rtrue),
     useBandInput = TRUE
   )
   testthat::expect_equal(length(out), 3)
})
