library(argo)
library(deSolve)

#### all possible data files

subdirs <- c("../results/for_paper/7param//variablephi-notemper", 
             "../results/for_paper/7param//variablephi-temper-warmstart",
             "../results/for_paper/7param//variablephi-temper-warmstart")
all_files <- lapply(subdirs, list.files)
all_files <- lapply(all_files, function(x) x[grep(".*log-([0-9]+)-7param.*", x)])
all_seeds <- lapply(all_files, function(x) gsub(".*log-([0-9]+)-7param.*", "\\1", x))
common_seeds <- all_seeds[[1]]
for (i in 2:length(all_seeds)){
  common_seeds <- intersect(common_seeds, all_seeds[[i]])  
}

common_seeds <- unique(all_seeds[[1]])

envhes1log <- new.env()


args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) == 0){
  rda_it = 1
}else{
  rda_it = args
}

##### Read incomplete dataset
each_seed <- common_seeds[rda_it]
dir.create("../results/for_paper/7param/variablephi-cool-0p33-warmstart", showWarnings = FALSE)

if (file.exists(paste0("../results/for_paper/7param/variablephi-cool-0p33-warmstart/Hes1-log-",
                       each_seed,"-7param-variablephi-cool-0p33-warmstart.rda"))){
  quit(save = "no")
}

load(paste0("../results/for_paper/7param/variablephi-notemper/Hes1-log-",each_seed,"-7param-variablephi-notemper.rda"))

# cooling up with (1, 1, 1/3) temperature ------------------------------------
#' there are 99 sampled X's (3 x 33) and only 33 actual observations.  
#' Using the logic about cooling, the default temperatures for Hes1 should be (1, 1, 1/3).
config$priorTemperature <- 1
priorTemperatureObs <- 1/3

xInit <- apply(gpode$xsampled, 2:3, mean)
thetaInit <- colMeans(gpode$theta)
phiNoTemperOptimized <- phiUsed

samplesCpp <- gpds:::solveGpdsRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = hes1logmodel,
  tvecFull = xsim$time,
  sigmaExogenous = pram.true$sigma,
  phiExogenous = phiNoTemperOptimized,
  xInitExogenous = xInit,
  thetaInitExogenous = thetaInit,
  muExogenous = matrix(numeric(0)),
  dotmuExogenous = matrix(numeric(0)),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = priorTemperatureObs,
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

phiUsed <- samplesCpp$phi
samplesCpp <- samplesCpp$llikxthetasigmaSamples

samplesCpp <- samplesCpp[,,1]
out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim[,-1]))
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(hes1logmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim[,-1]))

matplot(xsim$time, xCpp, type="l", add=TRUE, lty=2)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(hes1logmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::hes1logmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::hes1logmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- "../results/cpp/7param/variablephi-temper-warmstart/"
system(paste("mkdir -p", outDir))

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(phiUsed[,j], 2), collapse = "; ")
}

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-7param-variablephi-temper-warmstart.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save.image(paste0(outDir, config$modelName,"-",config$seed,"-7param-variablephi-temper-warmstart.rda"))
