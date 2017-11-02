#### first run without mean ----------------------------------------------------
library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123, #(as.integer(Sys.time())*132455+sample(1e6,1))%%1e7,
  npostplot = 5,
  loglikflag = "withmean",
  bandsize = 20,
  hmcSteps = 100,
  n.iter = 1e4,
  burninRatio = 0.1,
  stepSizeFactor = 0.1,
  bandwidthFactor = 0.5
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
             data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, config$kernel)
fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), 
                                  r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), 
                                            r.nobs, config$kernel)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)

factorList <- seq(0, 2, 0.01)[-1]
fn.sim.obs <- data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2])
fn.sim.obs.minusMu <- t(t(fn.sim.obs) - colMeans(fn.sim.obs))

lliksurface <- sapply(factorList, function(f){
  mypar <- marlikmap$par
  mypar[2] <- mypar[2]*f
  mypar[4] <- mypar[4]*f
  phisigllikC(mypar , fn.sim.obs.minusMu, 
              r.nobs, config$kernel)$value  
})
plot(factorList, lliksurface, type="l")
abline(v=0.5, col=2)

# marlikmap$par[c(1,3)] <- marlikmap$par[c(1,3)]*3^2
# marlikmap$par[c(2,4)] <- marlikmap$par[c(2,4)]*config$bandwidthFactor
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
cursigma <- marlikmap$par[5]
curCovV$mu <- rep(colMeans(fn.sim.obs)[1], length(curCovV$mu))
curCovR$mu <- rep(colMeans(fn.sim.obs)[2], length(curCovR$mu))

dotmu <- fODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
curCovV$dotmu <- rep(0, length(curCovV$mu))
curCovR$dotmu <- rep(0, length(curCovR$mu))

nall <- nrow(fn.sim)
burnin <- as.integer(config$n.iter*config$burninRatio)

xInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag)
set.seed(Sys.time())
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
gpds:::plotPostSamples(paste0("~/Desktop/",
                              config$loglikflag,"-experiment-",config$kernel,"-",config$seed,".pdf"), 
                       fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))
