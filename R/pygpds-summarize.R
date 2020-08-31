# Summarize the results fn.py ----
library(gpds)

config <- list()
config$modelName <- "FN"
config$noise <- rep(0.2, 2)

seeds <- list.files("~/Workspace/DynamicSys/dynamic-systems/pygpds/")  ## get the list of seeds ran
seeds <- seeds[grep(".*fn_inferred_trajectory_seed([0-9]+).*", seeds)]
seeds <- unique(gsub(".*fn_inferred_trajectory_seed([0-9]+).*", "\\1", seeds))

fnmodel <- list(
  fOde=gpds:::fODE,
  fOdeDx=gpds:::fnmodelDx,
  fOdeDtheta=gpds:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
}


theta_true <- c(0.2, 0.2, 3)
times <- seq(0, 20, length=241)
xtrue <- deSolve::ode(y = c(-1, 1), times = times, func = modelODE, parms = theta_true)
xtrue <- data.frame(xtrue)
rowId <- sapply(seq(0, 20, 0.5), function(x) which(abs(x-times) < 1e-6))
xtrue.obs <- xtrue[rowId, ]


oursPostX <- list()
oursPostTheta <- list()
oursPostSigma <- list()
oursRMSE <- list()

for(i in seeds){
  inferred_trajectory <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/fn_inferred_trajectory_seed", i, ".csv"), header = FALSE, sep=" ")
  inferred_trajectory <- t(inferred_trajectory)
  inferred_theta <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/fn_inferred_theta_seed", i, ".csv"), header = FALSE)
  inferred_theta <- unlist(inferred_theta)
  inferred_sigma <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/fn_inferred_sigma_seed", i, ".csv"), header = FALSE)
  oursPostX[[as.character(i)]]<- inferred_trajectory
  oursPostTheta[[as.character(i)]] <- inferred_theta
  oursPostSigma[[as.character(i)]] <- inferred_sigma
  xdesolve <- deSolve::ode(y = inferred_trajectory[1,], times = times, func = modelODE, parms = unlist(inferred_theta))
  rmse <- sqrt(colMeans((xdesolve[rowId, -1] - xtrue.obs[, -1])^2))
  oursRMSE[[as.character(i)]] <- rmse
}
oursRMSE <- do.call(rbind, oursRMSE)
oursPostTheta <- do.call(rbind, oursPostTheta)
colMeans(oursRMSE)
colMeans(oursPostTheta)


# Summarize the results of hes1.py
library(gpds)

config <- list()
config$modelName <- "Hes1log"
config$noise <- rep(0.2, 2)

seeds <- list.files("~/Workspace/DynamicSys/dynamic-systems/pygpds/")  ## get the list of seeds ran
seeds <- seeds[grep(".*hes1log_inferred_trajectory_seed([0-9]+).*", seeds)]
seeds <- unique(gsub(".*hes1log_inferred_trajectory_seed([0-9]+).*", "\\1", seeds))

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1logmodelODE(parameters, t(state))))
}

theta_true <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
times <- seq(0, 60*4, by = 0.01)
xtrue <- deSolve::ode(y = log(c(1.438575, 2.037488, 17.90385)), times = times, func = modelODE, parms = theta_true)
xtrue <- data.frame(xtrue)
rowId <- sapply(seq(0, 240, length=33), function(x) which(abs(x-times) < 1e-6))
xtrue.obs <- xtrue[rowId, ]

oursPostX <- list()
oursPostTheta <- list()
oursPostSigma <- list()
oursRMSE <- list()

for(i in seeds){
  inferred_trajectory <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/hes1log_inferred_trajectory_seed", i, ".csv"), header = FALSE, sep=" ")
  inferred_trajectory <- t(inferred_trajectory)
  inferred_theta <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/hes1log_inferred_theta_seed", i, ".csv"), header = FALSE)
  inferred_theta <- unlist(inferred_theta)
  inferred_sigma <- read.csv(paste0("~/Workspace/DynamicSys/dynamic-systems/pygpds/hes1log_inferred_sigma_seed", i, ".csv"), header = FALSE)
  oursPostX[[as.character(i)]]<- inferred_trajectory
  oursPostTheta[[as.character(i)]] <- inferred_theta
  oursPostSigma[[as.character(i)]] <- inferred_sigma
  xdesolve <- deSolve::ode(y = inferred_trajectory[1,], times = times, func = modelODE, parms = unlist(inferred_theta))
  rmse <- sqrt(colMeans((xdesolve[rowId, -1] - xtrue.obs[, -1])^2))
  oursRMSE[[as.character(i)]] <- rmse
}
oursRMSE <- do.call(rbind, oursRMSE)
oursPostTheta <- do.call(rbind, oursPostTheta)
colMeans(oursRMSE)
colMeans(oursPostTheta)


# Summarize the results of hes1-log-clean.R
library(gpds)

config <- list()
config$modelName <- "Hes1log"
config$noise <- rep(0.2, 2)

seeds <- list.files("~/Workspace/DynamicSys/results/for_paper/hes1log/variablephi-heating")  ## get the list of seeds ran
seeds <- seeds[grep(".*Hes1-log-([0-9]+)-hes1log_inferred_theta.csv", seeds)]
seeds <- unique(gsub(".*Hes1-log-([0-9]+)-hes1log_inferred_theta.csv", "\\1", seeds))

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1logmodelODE(parameters, t(state))))
}

theta_true <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
times <- seq(0, 60*4, by = 0.01)
xtrue <- deSolve::ode(y = log(c(1.438575, 2.037488, 17.90385)), times = times, func = modelODE, parms = theta_true)
xtrue <- data.frame(xtrue)
rowId <- sapply(seq(0, 240, length=33), function(x) which(abs(x-times) < 1e-6))
xtrue.obs <- xtrue[rowId, ]

oursPostX <- list()
oursPostTheta <- list()
oursRMSE <- list()

for(i in seeds){
  inferred_trajectory <- read.csv(paste0("~/Workspace/DynamicSys/results/for_paper/hes1log/variablephi-heating/Hes1-log-", i, "-hes1log_inferred_trajectory.csv"))
  inferred_trajectory <- data.matrix(inferred_trajectory[,-1])
  inferred_theta <- read.csv(paste0("~/Workspace/DynamicSys/results/for_paper/hes1log/variablephi-heating/Hes1-log-", i, "-hes1log_inferred_theta.csv"))
  inferred_theta <- unlist(inferred_theta[,-1])
  oursPostX[[as.character(i)]]<- inferred_trajectory
  oursPostTheta[[as.character(i)]] <- inferred_theta
  xdesolve <- deSolve::ode(y = inferred_trajectory[1,], times = times, func = modelODE, parms = unlist(inferred_theta))
  rmse <- sqrt(colMeans((xdesolve[rowId, -1] - xtrue.obs[, -1])^2))
  oursRMSE[[as.character(i)]] <- rmse
}
oursRMSE <- do.call(rbind, oursRMSE)
oursPostTheta <- do.call(rbind, oursPostTheta)
colMeans(oursRMSE)
colMeans(oursPostTheta)
