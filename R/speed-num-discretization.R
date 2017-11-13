library(testthat)
library(gpds)

# c++ profiling ---------------------------------------------------------------
config <- list(
  nobs = 41,
  noise = 0.1,
  seed = 123,
  bandsize = 20,
  kernel = "matern",
  filllevel = 4
)

VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=config$noise
)
fn.true <- VRtrue   # number of reference points is now flexible
fn.true$time <- seq(0,20,0.05)
fn.true.Vfunc <- approxfun(fn.true$time, fn.true$Vtrue)
fn.true.Rfunc <- approxfun(fn.true$time, fn.true$Rtrue)

fn.sim <- fn.true

set.seed(config$seed)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=config$noise)
fn.sim.obs <- fn.sim[seq(1,nrow(fn.sim), length=config$nobs),]
fn.sim <- insertNaN(fn.sim.obs, config$filllevel)

tvec.full <- fn.sim$time

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

curCovV$mu <- fn.true.Vfunc(fn.sim$time)  # pretend these are the means
curCovR$mu <- fn.true.Rfunc(fn.sim$time)

dotmu <- fODE(pram.true$abc, cbind(curCovV$mu, curCovR$mu)) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])

fn.sim.all <- fn.sim
lapseTime.all <- list()
n.all <- seq(17, 641, by=16)

for(nthis in n.all){
  fn.sim  <- head(fn.sim.all, nthis)
  
  tvec.full <- fn.sim$time
  
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
  
  curCovV$mu <- fn.true.Vfunc(fn.sim$time)  # pretend these are the means
  curCovR$mu <- fn.true.Rfunc(fn.sim$time)
  
  dotmu <- fODE(pram.true$abc, cbind(curCovV$mu, curCovR$mu)) # pretend these are the means for derivatives
  curCovV$dotmu <- as.vector(dotmu[,1])  
  curCovR$dotmu <- as.vector(dotmu[,2])
  
  lapseTime <- parallel::mclapply(1:16, function(dummy){
    speedRatio <- speedbenchmarkXthetallik( data.matrix(fn.sim[,1:2]),
                                            curCovV,
                                            curCovR,
                                            cursigma,
                                            c(sin(fn.sim$time), cos(fn.sim$time), 0.2, 0.2, 3),
                                            nrep = 1e3)
    rownames(speedRatio) <- c("xthetallik_rescaled", "xthetallikBandApproxHardCode",
                              "xthetallikHardCode", "xthetallik",
                              "xthetallik_withmu", "xthetallik_withmu2",
                              "xthetallikBandApprox", "xthetallikWithmuBand")
    speedRatio
  }, mc.cores=8)
  
  lapseTime <- do.call(cbind, lapseTime)
  
  lapseTime.all[[as.character(nthis)]] <- lapseTime
}
lapseTime.array <- sapply(lapseTime.all, identity, simplify = "array")


dimnames(lapseTime.array)
pdf("speed with num of discretization.pdf")
for(funname in dimnames(lapseTime.array)[[1]]){
  plot(as.numeric(dimnames(lapseTime.array)[[3]]), 
       colMeans(lapseTime.array[funname,,])/1e9, main = funname,
       xlab = "num. of discretization", ylab = "avg time in seconds for 1000 likelihood evals")
}
dev.off()


# R profiling -----------------------------------------------------------------
lapseTimeR.all <- list()
for(nthis in n.all){
  fn.sim  <- head(fn.sim.all, nthis)
  
  tvec.full <- fn.sim$time
  
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
  
  curCovV$mu <- fn.true.Vfunc(fn.sim$time)  # pretend these are the means
  curCovR$mu <- fn.true.Rfunc(fn.sim$time)
  
  dotmu <- fODE(pram.true$abc, cbind(curCovV$mu, curCovR$mu)) # pretend these are the means for derivatives
  curCovV$dotmu <- as.vector(dotmu[,1])  
  curCovR$dotmu <- as.vector(dotmu[,2])
  
  stepSize <- rep(0.00035, 2*nthis+3)*0.00001
  
  xthetaValues <- c(sin(fn.sim$time), cos(fn.sim$time), 0.2, 0.2, 3)
  lapseTime <- microbenchmark::microbenchmark(
    x <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xthetaValues, stepSize, 1e2, F, loglikflag = "withmeanBand")  
  )
  lapseTimeR.all[[as.character(nthis)]] <- lapseTime$time
}
lapseTimeR.all <- do.call(rbind, lapseTimeR.all)

plot(as.numeric(rownames(lapseTimeR.all)), 
     rowMeans(lapseTimeR.all)/1e9, main = "withmeanBand",
     xlab = "num. of discretization", ylab = "avg time in seconds for 100 leapfrog steps")
