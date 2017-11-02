#### second run with true mean -------------------------------------------------
library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e6,1))%%1e9,
  npostplot = 5,
  loglikflag = "withmean",
  bandsize = 20,
  hmcSteps = 200,
  n.iter = 1e4,
  burninRatio = 0.1,
  stepSizeFactor = 1,
  refitPhiSigma_withmu = FALSE,
  startAtTruth = TRUE,
  useTrueMu = TRUE
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

fn <- function(par) -phisigllikC( par, yobsGPsmooth, r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, yobsGPsmooth, r.nobs, config$kernel)$grad)
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
stepLow.traj <- xth.formal <- matrix(NA, n.iter, numparam)

xth.formal[1,] <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)  


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
  
  if(t %% 100 ==0) show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nall*2+1):(nall*2+3)]))
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
  "../results/2017-11-02/withmean_repSample_pluginTruth/",
  config$loglikflag,"-experiment-",config$kernel,"-",config$startAtTruth, "-", config$useTrueMu,".pdf"), 
  fn.true, fn.sim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  "../results/2017-11-02/withmean_repSample_pluginTruth/",
  config$loglikflag,"-experiment-",config$kernel,"-",config$startAtTruth, "-", config$useTrueMu,".rda"))
