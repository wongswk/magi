#### first run with flat mean --------------------------------------------------
library(gpds)
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = 0.1,
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 5,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 1000,
    n.iter = 300,
    burninRatio = 0.1,
    stepSizeFactor = 0.1,
    filllevel = 4
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
outDir <- with(config, paste0("~/Workspace/DynamicSys/results/", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

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

fn.sim <- insertNaN(fn.sim.obs,config$filllevel)

tvec.full <- fn.sim$time
tvec.nobs <- fn.sim.obs$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, config$noise), 
             data.matrix(fn.sim.obs[,1:2]), r.nobs, config$kernel)
fn <- function(par) -phisigllikC( par, data.matrix(fn.sim.obs[,1:2]), 
                                  r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim.obs[,1:2]), 
                                            r.nobs, config$kernel)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
cursigma <- marlikmap$par[5]

curCovV$mu[] <- mean(fn.sim.obs$Vtrue)
curCovR$mu[] <- mean(fn.sim.obs$Rtrue)

nall <- nrow(fn.sim)
burnin <- as.integer(config$n.iter*config$burninRatio)

xInit <- c(fn.true.Vfunc(fn.sim$time), fn.true.Rfunc(fn.sim$time), pram.true$abc)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

vInit <- getMeanCurve(fn.sim.obs$time, fn.sim.obs$Vtrue, fn.sim$time, 
                      t(marlikmap$par[1:2]), cursigma, config$kernel)
rInit <- getMeanCurve(fn.sim.obs$time, fn.sim.obs$Rtrue, fn.sim$time, 
                      t(marlikmap$par[3:4]), cursigma, config$kernel)
# TODO: add back the pilot idea to solve initial moving to target area problem
# or use tempered likelihood instead
# xInit <- c(vInit, rInit, rep(1,3))

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag)
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)

xth.formal <- chainSamplesOut[[1]]
lliklist <- chainSamplesOut[[2]]

n.iter <- config$n.iter

gpode <- list(abc=xth.formal[-(1:burnin), (nall*2+1):(nall*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=xth.formal[-(1:burnin), (nall+1):(nall*2)],
              vtrue=xth.formal[-(1:burnin), 1:nall],
              lp__=lliklist[-(1:burnin)],
              lglik=lliklist[-(1:burnin)])
gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-phase1.pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-phase1-",config$kernel,"-",config$seed,".rds"))


muV <- colMeans(gpode$vtrue)
muR <- colMeans(gpode$rtrue)
startTheta <- colMeans(gpode$abc)
dotmuV <- rowMeans(gpode$fode[,1,])
dotmuR <- rowMeans(gpode$fode[,2,])

#### second run with 1st phase mean -------------------------------------------

curCovV$mu <- muV
curCovR$mu <- muR

curCovV$dotmu <- dotmuV
curCovR$dotmu <- dotmuR

cursigma <- mean((data.matrix(fn.sim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(fn.sim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- config$n.iter

xInit <- c(muV, muR, startTheta)
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag)
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)


gpode <- list(abc=chainSamplesOut$xth[-(1:burnin), (nall*2+1):(nall*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=chainSamplesOut$xth[-(1:burnin), (nall+1):(nall*2)],
              vtrue=chainSamplesOut$xth[-(1:burnin), 1:nall],
              lp__=chainSamplesOut$lliklist[-(1:burnin)],
              lglik=chainSamplesOut$lliklist[-(1:burnin)])


gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-phase2.pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-phase2-",config$kernel,"-",config$seed,".rds"))

#### last run with true mean --------------------------------------------------

curCovV$mu <- fn.true.Vfunc(fn.sim$time)  # pretend these are the means
curCovR$mu <- fn.true.Rfunc(fn.sim$time)

dotmu <- fODE(pram.true$abc, cbind(curCovV$mu, curCovR$mu)) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])

cursigma <- mean((data.matrix(fn.sim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(fn.sim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- config$n.iter

xInit <- c(curCovV$mu, curCovR$mu, pram.true$abc)  
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag)
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)


gpode <- list(abc=chainSamplesOut$xth[-(1:burnin), (nall*2+1):(nall*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=chainSamplesOut$xth[-(1:burnin), (nall+1):(nall*2)],
              vtrue=chainSamplesOut$xth[-(1:burnin), 1:nall],
              lp__=chainSamplesOut$lliklist[-(1:burnin)],
              lglik=chainSamplesOut$lliklist[-(1:burnin)])

gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-trueMu.pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-trueMu-",config$kernel,"-",config$seed,".rds"))
