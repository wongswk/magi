#### first run with flat mean --------------------------------------------------
library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "matern",
  seed = (as.integer(Sys.time())*104729+sample(1e6,1))%%1e9,
  npostplot = 5,
  loglikflag = "withmean",
  bandsize = 20,
  hmcSteps = 200,
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

xInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

vInit <- getMeanCurve(fn.sim.obs$time, fn.sim.obs$Vtrue, fn.true$time, 
                      t(marlikmap$par[1:2]), cursigma, config$kernel)
rInit <- getMeanCurve(fn.sim.obs$time, fn.sim.obs$Rtrue, fn.true$time, 
                      t(marlikmap$par[3:4]), cursigma, config$kernel)
# TODO: add back the pilot idea to solve initial moving to target area problem
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

fn.sim$time <- fn.sim.all$time
gpds:::plotPostSamples(paste0(
  "~/Workspace/DynamicSys/results/2017-11-03/withmean_repSample_pluginPostMean/",
  config$loglikflag,"-experiment-",config$kernel,"-",config$seed,".pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

colMeans(gpode$abc)

muV <- colMeans(gpode$vtrue)
muR <- colMeans(gpode$rtrue)
startTheta <- colMeans(gpode$abc)
dotmuV <- rowMeans(gpode$fode[,1,])
dotmuR <- rowMeans(gpode$fode[,2,])

#### second run with mean ------------------------------------------------------
library(gpds)

config$loglikflag <- "withmean"
config$startAtTruth <- FALSE
config$useTrueMu <- FALSE

curCovV$mu <- muV
curCovR$mu <- muR

curCovV$dotmu <- dotmuV
curCovR$dotmu <- dotmuR

cursigma <- mean((data.matrix(fn.sim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(fn.sim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- config$n.iter
stepLow.traj <- xth.formal <- matrix(NA, n.iter, numparam)

xth.formal[1,] <- c(muV, muR, startTheta)

lliklist <- accepts <- c()
accepts[1] <- 1

burnin <- as.integer(n.iter*config$burninRatio)
stepLow <- rep(0.00035, 2*nall+3)*config$stepSizeFactor
timenow <- Sys.time()
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, config$hmcSteps, F, loglikflag = config$loglikflag)
  if(any(is.na(foo$final))) stop()
  xth.formal[t,] <- foo$final
  accepts[t] <- foo$acc
  stepLow.traj[t,] <- stepLow
  
  if (t < burnin & t > 10) {
    if (mean(tail(accepts[1:t],100)) > 0.9) {
      stepLow <- stepLow * 1.005
    } else if (mean(tail(accepts[1:t],100)) < 0.6) {
      stepLow <- stepLow * .995
    }
    xthsd <- apply(tail(xth.formal[1:t,],n.iter*config$burninRatio*0.1), 2, sd)
    if(mean(xthsd)>0) stepLow <- 0.01*xthsd/mean(xthsd)*mean(stepLow) + 0.99*stepLow
  }
  lliklist[t] <- foo$lpr
  
  if(t %% 100 ==0) methods::show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nall*2+1):(nall*2+3)]))
}

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

fn.sim$time <- fn.sim.all$time

gpds:::plotPostSamples(paste0(
  "~/Workspace/DynamicSys/results/2017-11-03/withmean_repSample_pluginPostMean/",
  config$loglikflag,"-experiment-",config$kernel,"-",config$seed,".pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  "~/Workspace/DynamicSys/results/2017-11-03/withmean_repSample_pluginPostMean/",
  config$loglikflag,"-experiment-",config$kernel,"-",config$seed,".rds"))
