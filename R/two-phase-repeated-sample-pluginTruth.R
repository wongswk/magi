#' deprecated on 2017-11-14, DO NOT RUN
#### second run with true mean -------------------------------------------------
library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e6,1))%%1e9,
  npostplot = 5,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 1000,
  n.iter = 1e4,
  burninRatio = 0.1,
  stepSizeFactor = 0.1,
  startAtTruth = TRUE,
  useTrueMu = TRUE,
  filllevel = 4
)
outDir <- "~/Workspace/DynamicSys/results/2017-11-08/withmean_repSample_pluginTruth/"
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
burnin <- as.integer(n.iter*config$burninRatio)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

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
  config$loglikflag,"-trueMu-",config$kernel,"-",config$seed,".pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-trueMu-",config$kernel,"-",config$seed,".rds"))
