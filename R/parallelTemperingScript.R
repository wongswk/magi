library(gpds)
config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123,
  npostplot = 5,
  loglikflag = "band",
  bandsize = 20,
  hmcSteps = 100, # FIXME put this to 1000 program will not finish
  n.iter = 1e4,
  burninRatio = 0.1,
  stepSizeFactor = 1
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

marlikmap <- list(par=c(2.314334, 1.346233, 0.622316, 2.451729, 0.084745))

cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)

curCovV$mu <- as.vector(fn.true[,1])  # pretend these are the means
curCovR$mu <- as.vector(fn.true[,2])

dotmu <- fODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])


nall <- nrow(fn.sim)
numparam <- nall*2+3

burnin <- as.integer(config$n.iter*config$burninRatio)
stepLow <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

out <- parallel_temper_hmc_xtheta(
  data.matrix(fn.sim[,1:2]),
  curCovV,
  curCovR,
  cursigma,
  3^(seq(0,1,length=7)),
  0.1,
  c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc),
  stepLow,
  config$hmcSteps,
  config$n.iter)

plot(out[405,1,], type="l")
plot(gpode$abc[,2], type="l")

hist(out[405,1,])
hist(out[406,1,])

hist(gpode$abc[,1], breaks = 40, col=rgb(1,0,0,0.5), probability = T)
hist(out[404,1,], breaks = 40, col=rgb(0,1,0,0.5), add=T, probability = T)

hist(gpode$abc[,2], breaks = 40, col=rgb(1,0,0,0.5), probability = T)
hist(out[405,1,], breaks = 40, col=rgb(0,1,0,0.5), add=T, probability = T)

hist(gpode$abc[,3], breaks = 40, col=rgb(1,0,0,0.5), probability = T)
hist(out[406,1,], breaks = 40, col=rgb(0,1,0,0.5), add=T, probability = T)
