# This script runs one replication of the HIV time-dependent model with randomized values of the hyperparameters phi

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) == 0){
  seed = 0
}else{
  seed = args
}

library(magi)

# HIV time-dependent model

hivtdmodelODE <- function(theta, x, tvec) {
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]
  
  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))
  
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] = lambda - rho * TU - eta * TU * V
  result[, 2] = eta * TU * V - delta * TI
  result[, 3] = N * delta * TI - c * V
  
  result
}

hivtdmodelDx <- function(theta, x, tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  lambda <- theta[1]
  rho <- theta[2]
  delta <- theta[3]
  N <- theta[4]
  c <- theta[5]
  
  eta <- 9e-5 * (1 - 0.9 * cos(pi * tvec / 1000))
  
  resultDx[, , 1] = cbind(-rho - eta * V, 0, -eta * TU)
  resultDx[, , 2] = cbind(eta * V, -delta, eta * TU)
  resultDx[, , 3] = cbind(rep(0, nrow(x)), N * delta, -c)
  
  resultDx
}

hivtdmodelDtheta <- function(theta, x, tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  TU <- x[, 1]
  TI <- x[, 2]
  V <- x[, 3]
  
  delta <- theta[3]
  N <- theta[4]
  
  resultDtheta[, , 1] = cbind(1, -TU, 0, 0, 0)
  resultDtheta[, , 2] = cbind(0, 0, -TI, 0, 0)
  resultDtheta[, , 3] = cbind(0, 0, N * TI, delta * TI, -V)
  
  resultDtheta
}

hivtdmodel <- list(
  fOde = hivtdmodelODE,
  fOdeDx = hivtdmodelDx,
  fOdeDtheta = hivtdmodelDtheta,
  thetaLowerBound = rep(0, 5),
  thetaUpperBound = rep(Inf, 5)
)

param.true <- list(
  theta = c(36, 0.108, 0.5, 1000, 3),
  x0 = c(600, 30, 1e5),
  sigma = c(sqrt(10), sqrt(10), 10),
  times = seq(0, 20, 0.2)
)

set.seed(12321)
modelODE <- function(tvec, state, parameters) {
  list(as.vector(hivtdmodelODE(parameters, t(state), tvec)))
}

xtrue <- deSolve::ode(y = param.true$x0, times = param.true$times,
                      func = modelODE, parms = param.true$theta)
y <- data.frame(xtrue)
for(j in 1:(ncol(y) - 1)){
  y[, 1+j] <- y[, 1+j] + rnorm(nrow(y), sd = param.true$sigma[j])
}

compnames <- c("TU", "TI", "V")
complabels <- c("Concentration", "Concentration", "Load")

phiEst <- matrix(0, nrow = 2, ncol = ncol(y) - 1)
sigmaInit <- rep(0, ncol(y) - 1)
for (j in 1:(ncol(y) - 1)){
  hyperparam <- gpsmoothing(y[, j+1], y[, "time"])
  phiEst[, j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
}

colnames(phiEst) <- compnames
phiEst
sigmaInit

phiEst[, 3] <- c(1e7, 0.5)
sigmaInit[3] <- 100

y_I <- setDiscretization(y, level = 1)

if (seed > 0) {
  set.seed(seed)             
  for (j in 1:3) {
    phiEst[, j] <- runif(2, 2/3 * phiEst[, j], 3/2 * phiEst[, j])
  }                                  
}

HIVresult <- MagiSolver(y_I, hivtdmodel,
                        control = list(phi = phiEst, sigma = sigmaInit, verbose = TRUE))

save(HIVresult, y_I, y, file=paste0("random-phi-", seed, ".rda"))
