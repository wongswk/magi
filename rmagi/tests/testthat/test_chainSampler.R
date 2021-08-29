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

r.nobs <- abs(outer(tvec.nobs, t(tvec.nobs),'-')[,1,])
yobs1 <- data.matrix(fn.sim.obs[,1,drop=FALSE])
outputc <- magi:::gpsmooth(yobs1,
                           r.nobs,
                           config$kernel, sigmaExogenScalar = -1, FALSE)

