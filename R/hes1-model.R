#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
if(!exists("config")){
  config <- list(
    nobs = 51,
    noise = 0.1,
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 5,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 200,
    n.iter = 1000,
    burninRatio = 0.1,
    stepSizeFactor = 0.1,
    filllevel = 2,
    modelName = "hes1"
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise",noise,"-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  xinit = c(1, 1, 1),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 30, by = 0.001)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1modelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = c(1, 1, 1), times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- xtrue

set.seed(config$seed)
xsim[,-1] <- xsim[,-1]+rnorm(length(unlist(xsim[,-1])), sd=config$noise)
xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]

xsim <- insertNaN(xsim.obs,config$filllevel)

tvec.full <- xsim$time
tvec.nobs <- xsim.obs$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, 1, 1, config$noise), 
             data.matrix(xsim.obs[,-1]), r.nobs, config$kernel)
fn <- function(par) -phisigllikC( par, data.matrix(xsim.obs[, -1]), 
                                  r.nobs, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(xsim.obs[, -1]),
                                            r.nobs, config$kernel)$grad)
marlikmap <- optim(rep(1, 2*ncol(xsim.obs)-1), fn, gr, method="L-BFGS-B", lower = 0.0001)

cursigma <- tail(marlikmap$par, 1)

curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
  covEach <- calCov(marlikmap$par[(2*j-1):(2*j)], r, signr, bandsize=config$bandsize, 
                    kerneltype=config$kernel)
  covEach$mu[] <- mean(xsim.obs[,j+1])
  covEach
})

curCovV$mu[] <- mean(xsim.obs$Vtrue)
curCovR$mu[] <- mean(xsim.obs$Rtrue)

nall <- nrow(xsim)
burnin <- as.integer(config$n.iter*config$burninRatio)

xInit <- c(xtrue.Vfunc(xsim$time), xtrue.Rfunc(xsim$time), pram.true$abc)
stepLowInit <- rep(0.00035, 2*nall+3)*config$stepSizeFactor

vInit <- getMeanCurve(xsim.obs$time, xsim.obs$Vtrue, xsim$time, 
                      t(marlikmap$par[1:2]), cursigma, config$kernel)
rInit <- getMeanCurve(xsim.obs$time, xsim.obs$Rtrue, xsim$time, 
                      t(marlikmap$par[3:4]), cursigma, config$kernel)
# TODO: add back the pilot idea to solve initial moving to target area problem
# or use tempered likelihood instead
# xInit <- c(vInit, rInit, rep(1,3))

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,1:2]), list(curCovV, curCovR), cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature)
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

xtrue$dVtrue = with(c(xtrue,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
xtrue$dRtrue = with(c(xtrue,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-priorTempered.pdf"), 
  xtrue, xsim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTempered-",config$kernel,"-",config$seed,".rds"))


muV <- colMeans(gpode$vtrue)
muR <- colMeans(gpode$rtrue)
startTheta <- colMeans(gpode$abc)
dotmuV <- rowMeans(gpode$fode[,1,])
dotmuR <- rowMeans(gpode$fode[,2,])

# run with priorTempered phase 2 -------------------------------------------

curCovV$mu <- muV
curCovR$mu <- muR

curCovV$dotmu <- dotmuV
curCovR$dotmu <- dotmuR

cursigma <- mean((data.matrix(xsim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(xsim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(xsim)
numparam <- nall*2+3
n.iter <- config$n.iter

xInit <- c(muV, muR, startTheta)
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,1:2]), list(curCovV, curCovR), cursigma, 
               xthetaValues, stepSize, config$hmcSteps, F, loglikflag = config$loglikflag,
               priorTemperature = config$priorTemperature)
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

xtrue$dVtrue = with(c(xtrue,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
xtrue$dRtrue = with(c(xtrue,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-priorTemperedPhase2.pdf"), 
  xtrue, xsim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-priorTemperedPhase2-",config$kernel,"-",config$seed,".rds"))

#### last run with true mean --------------------------------------------------

curCovV$mu <- xtrue.Vfunc(xsim$time)  # pretend these are the means
curCovR$mu <- xtrue.Rfunc(xsim$time)

dotmu <- fODE(pram.true$abc, cbind(curCovV$mu, curCovR$mu)) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])

cursigma <- mean((data.matrix(xsim[,1:2])-cbind(curCovV$mu, curCovR$mu))[!is.nan(xsim[,1]),]^2)
cursigma <- sqrt(cursigma)

nall <- nrow(xsim)
numparam <- nall*2+3
n.iter <- config$n.iter

xInit <- c(curCovV$mu, curCovR$mu, pram.true$abc)  
stepLowInit <- chainSamplesOut$stepLow

singleSampler <- function(xthetaValues, stepSize) 
  xthetaSample(data.matrix(xsim[,1:2]), list(curCovV, curCovR), cursigma, 
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

xtrue$dVtrue = with(c(xtrue,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
xtrue$dRtrue = with(c(xtrue,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

gpds:::plotPostSamples(paste0(
  outDir,
  config$kernel,"-",config$seed,"-trueMu.pdf"), 
  xtrue, xsim, gpode, pram.true, config)

absCI <- apply(gpode$abc, 2, quantile, probs = c(0.025, 0.5, 0.975))
absCI <- rbind(absCI, mean=colMeans(gpode$abc))
absCI <- rbind(absCI, coverage = (absCI["2.5%",] < pram.true$abc &  pram.true$abc < absCI["97.5%",]))

saveRDS(absCI, paste0(
  outDir,
  config$loglikflag,"-trueMu-",config$kernel,"-",config$seed,".rds"))
