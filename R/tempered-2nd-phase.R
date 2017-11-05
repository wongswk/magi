load("/Users/shihaoyang/Workspace/DynamicSys/results/2017-10-29/withmean start mu/first_phase_inference.rda")
library(gpds)

config <- list(
  nobs = 41,
  noise = 0.1,
  kernel = "generalMatern",
  seed = 123,
  npostplot = 5,
  loglikflag = "withmean",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 10000,
  burninRatio = 0.1,
  stepSizeFactor = 0.005,
  refitPhiSigma_withmu = FALSE,
  startAtTruth = TRUE,
  useTrueMu = TRUE,
  temperatureT = 2000
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

fn <- function(par) -phisigllikC( par, yobsGPsmooth, 
                                  r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, yobsGPsmooth, 
                                            r.nobs, config$kernel)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
cursigma <- marlikmap$par[5]
if(config$useTrueMu){
  curCovV$mu <- as.vector(fn.true[,1])  # pretend these are the means
  curCovR$mu <- as.vector(fn.true[,2])
  
  dotmu <- fODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
  curCovV$dotmu <- as.vector(dotmu[,1])  
  curCovR$dotmu <- as.vector(dotmu[,2])
}else{
  curCovV$mu <- muV
  curCovR$mu <- muR
  
  curCovV$dotmu <- dotmuV
  curCovR$dotmu <- dotmuR
}


cursigma <- mean((data.matrix(fn.sim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(fn.sim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- config$n.iter
if(config$startAtTruth){
  xInit <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)  
}else{
  xInit <- c(muV, muR, startTheta)
}


config$temperatureT <- 1e4
config$stepSizeFactor <- 0.01

stepLowInit <- rep(0.01, 2*nall+3)*config$stepSizeFactor



singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag, 
               temperature = config$temperatureT)
chainSamplesOut <- chainSampler(config, xInit, singleSampler, stepLowInit, verbose=TRUE)

plot(chainSamplesOut$xth[,404], type="l")
hist(chainSamplesOut$xth[,404])

mean(chainSamplesOut$xth[-(1:with(config,n.iter*burninRatio)),404])

burnin <- with(config,n.iter*burninRatio)
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

fn.sim$time <- fn.sim.all$time
gpds:::plotPostSamples(paste0(config$loglikflag,"-tempered-2nd-phase-",config$kernel,".pdf"), 
                       fn.true, fn.sim, gpode, pram.true, config)


# investigate group move, and the temperature at which the system becomes white noise, 
# do a temperature truncation. 