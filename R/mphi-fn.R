# sample from $\| f - m X\|_I$ --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = TRUE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1e5),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# remove force $K$ to be identity matrix from 1 --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1e5),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# plug-in true curve as mean --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = TRUE,
  forsePhase8Mean = FALSE,
  priorTemperature = c(1, 1),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# plug-in biased phase8 inference curve as mean --------------------------------------------------------
rm(list=ls())
config <- list(
  nobs = 201,
  noise = c(0.5, 0.5),
  kernel = "generalMatern",
  mphiType = "zero",
  forceDiagKphi = FALSE,
  forseTrueMean = FALSE,
  forsePhase8Mean = TRUE,
  priorTemperature = c(1, 1),
  seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 2500,
  burninRatio = 0.10,
  stepSizeFactor = 1,
  filllevel = 0,
  modelName = "FN",
  startXAtTruth = TRUE,
  startThetaAtTruth = TRUE,
  startSigmaAtTruth = TRUE,
  useGPmean = TRUE,
  phase2 = FALSE,
  temperPrior = TRUE,
  phase3 = FALSE,
  max.epoch = 10,
  epoch_method = c("mean", "median", "deSolve", "f_x_bar")[4]
)
source("R/fn-model.R")

# STAN sampling using finite difference ----------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())
init <- list(
  abc = pram.true$theta,
  vtrue = xtrue.atsim[,1],
  rtrue = xtrue.atsim[,2]
)

stanConfig <- list(
  sigma_obs=0.5,
  sigma_xdot=0.1
)
config$stanConfig <- stanConfig

gpsmooth <- stan(file="stan/m-finiteDifference.stan",
                 data=list(
                   N=nrow(xsim),
                   vobs=xsim$X1,
                   robs=xsim$X2,
                   time=xsim$time,
                   sigma_obs=stanConfig$sigma_obs,
                   sigma_xdot=stanConfig$sigma_xdot
                 ),
                 iter=10000, chains=7, init=rep(list(init), 7), warmup = 1000)

gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

stanode <- list(theta=gpsmooth_ss$abc,
                xsampled=abind::abind(gpsmooth_ss$vtrue, gpsmooth_ss$rtrue, along = 3),
                lglik=gpsmooth_ss$lp__,
                sigma = pram.true$sigma)
stanode$fode <- sapply(1:length(stanode$lglik), function(t) 
  with(stanode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
stanode$fode <- aperm(stanode$fode, c(3,1,2))


gpds:::plotPostSamplesFlex(
  paste0("m-finisteDifference-stan.pdf"), 
  xtrue, dotxtrue, xsim, stanode, pram.true, config, odemodel = odemodel)

save.image(paste0(filename, "-stan.rda"))

# likelihood move away from mode ------------------------------------------------
gpsmooth <- stan(file="stan/m-finiteDifference.stan",
                 data=list(
                   n_discret=nrow(xsim),
                   n_obs=nrow(xsim),
                   vobs=xsim$X1,
                   robs=xsim$X2,
                   time=xsim$time,
                   obs_index=1:nrow(xsim),
                   sigma_obs=stanConfig$sigma_obs,
                   sigma_xdot=stanConfig$sigma_xdot
                 ),
                 iter=200, chains=1, init=rep(list(init), 1), warmup = 20)

gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

stanode <- list(theta=gpsmooth_ss$abc,
                xsampled=abind::abind(gpsmooth_ss$vtrue, gpsmooth_ss$rtrue, along = 3),
                lglik=gpsmooth_ss$lp__,
                sigma = pram.true$sigma)
stanode$fode <- sapply(1:length(stanode$lglik), function(t) 
  with(stanode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
stanode$fode <- aperm(stanode$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0("m-finisteDifference-likelihoodmove-stan.pdf"), 
  xtrue, dotxtrue, xsim, stanode, pram.true, config)


