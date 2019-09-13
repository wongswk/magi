#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 15,
    noise = rep(0.01, 5), # 0.001 = low noise, 0.01 = high noise
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 40,
    hmcSteps = 100,
    n.iter = 20001,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    #filllevel = 3,
    linfillspace = 0.1,
    t.end = 100,
    modelName = "PTrans",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    linearizexInit = TRUE,
    useExoSigma = TRUE,
    max.epoch = 5
  )
}

#config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$ndis <- config$t.end / config$linfillspace + 1;
if(config$temperPrior){
  config$priorTemperature <- config$ndis / config$nobs  
}else{
  config$priorTemperature <- 1
}

if(config$loglikflag == "withmeanBand"){
  config$useMean = TRUE
  config$useBand = TRUE
}else if(config$loglikflag == "band"){
  config$useMean = FALSE
  config$useBand = TRUE
}else if(config$loglikflag == "withmean"){
  config$useMean = TRUE
  config$useBand = FALSE
}else if(config$loglikflag == "usual"){
  config$useMean = FALSE
  config$useBand = FALSE
}

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,100,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::ptransmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

#xsim <- insertNaN(xsim.obs,config$filllevel)
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}

if (config$useExoSigma) {
  exoSigma = rep(0.001, ncol(xsim)-1) #config$noise
} else {
  exoSigma = numeric(0)
}

if (config$linearizexInit) {
  exoxInit <- sapply(2:ncol(xsim.obs), function(j)
    approx(xsim.obs[, "time"], xsim.obs[, j], xsim[, "time"])$y)
} else {
  exoxInit <- matrix(nrow=0,ncol=0)
}


# cpp inference ----------------------------
ptransmodel <- list(
  fOde=gpds:::ptransmodelODE,
  fOdeDx=gpds:::ptransmodelDx,
  fOdeDtheta=gpds:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(100,6)
)

samplesCppAll <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = ptransmodel,
  tvecFull = xsim$time,
  sigmaExogenous = exoSigma,
  phiExogenous = matrix(nrow=0,ncol=0),
  xInitExogenous = exoxInit,
  muExogenous = matrix(nrow=0,ncol=0),
  dotmuExogenous = matrix(nrow=0,ncol=0),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = config$n.iter,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

for (epoch in 1:config$max.epoch) {

  samplesCpp <- samplesCppAll[,,epoch]
  out <- samplesCpp[-1,1,drop=FALSE]
  xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim[,-1]))
  thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(ptransmodel$thetaLowerBound)), 1]
  sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))
  
  matplot(xsim$time, xCpp, type="l", add=TRUE)
  
  llikId <- 1
  xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
  thetaId <- (max(xId)+1):(max(xId)+length(ptransmodel$thetaLowerBound))
  sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))
  
  
  burnin <- as.integer(config$n.iter*config$burninRatio)
  gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
                xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                               dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
                lglik=samplesCpp[llikId,-(1:burnin)],
                sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
  gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
    with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  gpode$fode <- aperm(gpode$fode, c(3,1,2))
  
  dotxtrue = gpds:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))
  
  odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)
  
  outDir <- "../results/cpp/"
  
  gpds:::plotPostSamplesFlex(
    paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "-epoch", epoch, ".pdf"), 
    xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

}
