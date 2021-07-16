library(magi)

outDir <- "../results/hiv-time-dependent/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

hivtdmodelODE <- function(theta,x,tvec) {
  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = lambda - rho * TU - eta * TU * V
  result[,2] = eta * TU * V - delta * TI
  result[,3] = N * delta * TI - c * V

  result
}

hivtdmodelDx <- function(theta,x,tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))

  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  resultDx[,1,1] = -rho - eta * V
  resultDx[,2,1] = 0
  resultDx[,3,1] = -eta * TU

  resultDx[,1,2] = eta * V
  resultDx[,2,2] = -delta
  resultDx[,3,2] = eta * TU

  resultDx[,1,3] = 0
  resultDx[,2,3] = N * delta
  resultDx[,3,3] = -c

  resultDx
}

hivtdmodelDtheta <- function(theta,x,tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))

  TU <- x[,1]
  TI <- x[,2]
  V <- x[,3]

  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]

  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))

  resultDtheta[,1,1] = 1
  resultDtheta[,2,1] = -TU
  resultDtheta[,3,1] = 0
  resultDtheta[,4,1] = 0
  resultDtheta[,5,1] = 0

  resultDtheta[,1,2] = 0
  resultDtheta[,2,2] = 0
  resultDtheta[,3,2] = -TI
  resultDtheta[,4,2] = 0
  resultDtheta[,5,2] = 0

  resultDtheta[,1,3] = 0
  resultDtheta[,2,3] = 0
  resultDtheta[,3,3] = N * TI
  resultDtheta[,4,3] = delta * TI
  resultDtheta[,5,3] = -V

  resultDtheta
}


# set up configuration if not already exist ------------------------------------
config <- list(
  nobs = 101,
  noise = c(sqrt(10), sqrt(10), 10),
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  niterHmc = 20000,
  filllevel = 1,
  t.end = 20,
  t.interval = 0.1,
  modelName = "hiv-time-dependent"
)


# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(36, 0.108, 0.5, 1000, 3), # lambda, rho, delta, N, c
  x0 = c(600, 30, 1e5), # TU, TI, V
  phi = cbind(c(190^2, 4.5), c(107^2, 3), c(2e4^2, 1)),
  sigma=config$noise
)

times <- seq(0,config$t.end,config$t.interval)

modelODE <- function(t, state, parameters) {
  list(as.vector(hivtdmodelODE(parameters, t(state), t)))
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

xsim <- setDiscretization(xsim.obs,config$filllevel)

hivtdmodel <- list(
  fOde=hivtdmodelODE,
  fOdeDx=hivtdmodelDx,
  fOdeDtheta=hivtdmodelDtheta,
  thetaLowerBound=c(0,0,0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf,Inf,Inf)
)

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

testDynamicalModel(hivtdmodelODE, hivtdmodelDx, hivtdmodelDtheta, "HIV time-dependent system", xInitExogenous, pram.true$theta, xsim$time)

phiExogenous <- matrix(0, nrow=2, ncol=ncol(xsim)-1)
sigmaInit <- rep(0, ncol(xsim)-1)
for (j in 1:(ncol(xsim)-1)){
  hyperparam <- gpsmoothing(xsim.obs[,j+1],
                            xsim.obs$time,
                            "generalMatern")
  phiExogenous[,j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
}


#' manually override estimated hyper-parameters for component 3 because 
#' GP smoothing gives bad result for rapidly decreasing curve
phiExogenous[2,3] <- 1
sigmaInit[3] <- 10

OursStartTime <- proc.time()[3]  # hivtdmodel is implemented in R only, and thus the running will be slow 

result <- magi::MagiSolver(xsim[,-1], hivtdmodel, xsim$time, control = list(xInit = xInitExogenous, niterHmc=config$niterHmc, stepSizeFactor = 0.06, phi=phiExogenous, sigma=sigmaInit))

OursTimeUsed <- proc.time()[3] - OursStartTime

gpode <- result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, hivtdmodelODE(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = hivtdmodelODE(pram.true$theta, data.matrix(xtrue[,-1]), xtrue$time)

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}

gpode$lglik <- gpode$lp
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-hivtd_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-hivtd_inferred_theta.csv"))
