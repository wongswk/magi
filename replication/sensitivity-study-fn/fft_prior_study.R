library(magi)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) > 0){
  filllevel <- args[1]
  seed <- scan("fn-seeds.txt")[args[2]]
  nobs_keep <- args[3]
  time_horizon <- as.numeric(args[4])
  n_interpolation <- as.numeric(args[5])
}else{
  filllevel <- 4
  seed <- 117015794
  nobs_keep <- 21
  time_horizon <- 20
  n_interpolation <- 81
}

outDir <- paste0("../results/fn-sparse-fill", filllevel, "-nobs", nobs_keep, "-ninterp", n_interpolation, "-timeend", time_horizon, "/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = nobs_keep,
    noise = c(0.2, 0.2),
    kernel = "generalMatern",
    seed = seed,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = filllevel,
    t.end = time_horizon,
    modelName = "FN",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    max.epoch = 1
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperatureObs <- 1

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
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::fnmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,config$t.end,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- insertNaN(xsim.obs,config$filllevel)

fnmodel <- list(
  fOde=magi:::fODE,
  fOdeDx=magi:::fnmodelDx,
  fOdeDtheta=magi:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

config$priorTemperature <- config$ndis / config$nobs  
config$priorTemperatureObs <- 1

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

OursStartTime <- proc.time()[3]

# set phi initialization with n_interpolation points
for(n_interpolation in c(21, 41, 81, 161)){
  cat("++++++ n_interpolation =", n_interpolation, "+++++++\n")
  id4phi <- seq(1, nrow(xInitExogenous), length.out = n_interpolation)
  samplesCpp <- magi:::solveMagiRcpp(
    yFull = xInitExogenous[id4phi,],
    odeModel = fnmodel,
    tvecFull = xsim$time[id4phi],
    sigmaExogenous = numeric(0),
    phiExogenous = matrix(nrow=0,ncol=0),
    xInitExogenous = xInitExogenous[id4phi,],
    thetaInitExogenous = matrix(nrow=0,ncol=0),
    muExogenous = matrix(nrow=0,ncol=0),
    dotmuExogenous = matrix(nrow=0,ncol=0),
    priorTemperatureLevel = config$priorTemperature,
    priorTemperatureDeriv = config$priorTemperature,
    priorTemperatureObs = config$priorTemperatureObs,
    kernel = config$kernel,
    nstepsHmc = config$hmcSteps,
    burninRatioHmc = config$burninRatio,
    niterHmc = 2,
    stepSizeFactorHmc = config$stepSizeFactor,
    nEpoch = config$max.epoch,
    bandSize = config$bandsize,
    useFrequencyBasedPrior = config$useFrequencyBasedPrior,
    useBand = config$useBand,
    useMean = config$useMean,
    useScalerSigma = config$useScalerSigma,
    useFixedSigma = config$useFixedSigma,
    verbose = TRUE)
}


for(n_interpolation in c(21, 41, 81, 161)){
  cat("++++++ n_interpolation =", n_interpolation, "+++++++\n")
  id4phi <- seq(1, nrow(xInitExogenous), length.out = n_interpolation)
  priorFactorV <- magi:::getFrequencyBasedPrior(xInitExogenous[id4phi,1], showplot = TRUE)
  priorFactorR <- magi:::getFrequencyBasedPrior(xInitExogenous[id4phi,2])
  print("priorFactorV")
  print(priorFactorV)
  print("priorFactorR")
  print(priorFactorR)
}

pdf("fft_problem_illustration.pdf", height = 7, width = 10)
for(n_interpolation in c(21, 41, 81, 161, 321, 641)){
  cat("++++++ n_interpolation =", n_interpolation, "+++++++\n")
  id4phi <- seq(1, nrow(xInitExogenous), length.out = n_interpolation)
  par(mar=c(5,3,11,3))
  plot(xsim$time[id4phi], xInitExogenous[id4phi,1], main=paste0("component V with n_interpolation = ", n_interpolation), xlab="time")
  axis(3, at=5*(0:4), labels = seq(1, length(id4phi), length.out = 5))
  mtext("index")
  
  x <- xInitExogenous[id4phi,1]
  z <- fft(x)
  zmod <- Mod(z)
  names(zmod) <- NULL
  plot(zmod, col=c(1,rep(2, (length(zmod)-1)/2),rep(1, (length(zmod)-1)/2)))
  title("modulus of fft", line=7)
  
  zmodEffective <- zmod[-1]
  zmodEffective <- zmodEffective[1:(length(zmodEffective)/2)]
  names(zmodEffective) <- 1:length(zmodEffective)
  
  upperQuarter = sort(zmodEffective)[ceiling(length(zmodEffective) * 0.75)]
  lowerQuarter = sort(zmodEffective)[floor(length(zmodEffective) * 0.25)]
  iqr = upperQuarter - lowerQuarter
  outliers = zmodEffective[zmodEffective > upperQuarter + 1.5 * iqr]
  if (length(outliers) == 0) {
    freq <- which.max(zmodEffective)
    abline(v=xsim$time[id4phi][freq], col=4)
  } else {
    outliers <- outliers[outliers > median(zmodEffective)]
    whichOutliers <- which(zmodEffective %in% outliers)
    abline(v=1+whichOutliers, col=4)
    freq <- max(whichOutliers)
  }
  meanFactor <- 0.5/freq
  legend("top", c("range of frequencies with large loadings", 
                  "modulus of effective frequency loading"), 
         lty=c(1, NA), pch=c(NA, 1), col=c(4, 2))
  msg <- paste0("WANT: frequency with largest loading = ", which.max(zmodEffective), 
                ", corresponding prior factor = ", round(0.5/which.max(zmodEffective), digits = 4))
  freq_wtd <- weighted.mean(whichOutliers, outliers^2)
  msg <- paste0(msg, "\nWANT: mod^2 weighted average frequency with large loadings = ", freq_wtd, ", corresponding prior factor = ", round(0.5/freq_wtd, digits = 4))
  freq_wtd <- weighted.mean(1:length(zmodEffective), zmodEffective^2)
  msg <- paste0(msg, "\nWANT: mod^2 weighted average frequency among all = ", freq_wtd, ", corresponding prior factor = ", round(0.5/freq_wtd, digits = 4))
  period_wtd <- weighted.mean(0.5/whichOutliers, outliers^2)
  msg <- paste0(msg, "\nWANT: mod^2 weighted average half periodicity (i.e. prior factor) with large loadings = ", period_wtd)
  period_wtd <- weighted.mean(0.5/(1:length(zmodEffective)), zmodEffective^2)
  msg <- paste0(msg, "\nWANT: mod^2 weighted average half periodicity (i.e. prior factor) among all = ", period_wtd)
  msg <- paste0(msg, "\nCURRENT: highest frequency with large loadings = ", freq, ", corresponding prior factor = ", round(0.5/freq, digits = 4))
  mtext(msg)
}
dev.off()

