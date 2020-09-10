library(magi)
library(testthat)

#when step size is too large, withmean gives NaN

context("x-theta sampler")
config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123,
  npostplot = 5,
  loglikflag = "band",
  bandsize = 20,
  hmcSteps = 100,
  n.iter = 1e3,
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



foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

marlikmap <- list(par = c(2.314, 1.346, 0.622, 2.452, 0.085))
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

xth <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)
rstep <- rep(1e10, length(xth))
foo <- xthetaSample(data.matrix(fn.sim[,1:2]), list(curCovV, curCovR), cursigma, 
                    xth, rstep, config$hmcSteps, F, loglikflag = "usual")
test_that("NA should not appear after xthetaSample", {
  expect_true(all(!is.na(foo$final)))
})

rstep <- rep(1e30, length(xth))
foo <- xthetaSample(data.matrix(fn.sim[,1:2]), list(curCovV, curCovR), cursigma, 
                    xth, rstep, config$hmcSteps, T, loglikflag = "usual")
test_that("NA should not appear after xthetaSample", {
  expect_true(all(!is.na(foo$final)))
})

rstep <- rep(0.0001, length(xth))
set.seed(234)
foo1 <- xthetaSample(data.matrix(fn.sim[,1:2]), list(curCovV, curCovR), cursigma, 
                    xth, rstep, config$hmcSteps, F, loglikflag = "usual")
set.seed(234)
foo2 <- xthetaSample(data.matrix(fn.sim[,1:2]), list(curCovV, curCovR), cursigma, 
                     xth, rstep, config$hmcSteps, T, loglikflag = "usual")
test_that("return trajectory or not does not affect result", {
  expect_equal(foo1$final, foo2$final)
})

