library("magi")

# Hes1 example

hes1modelODE <- function(theta, x, tvec) {
  P = x[, 1]
  M = x[, 2]
  H = x[, 3]
  
  PMHdt = array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * P * H + theta[2] * M - theta[3] * P
  PMHdt[, 2] = -theta[4] * M + theta[5] / (1 + P^2)
  PMHdt[, 3] = -theta[1] * P * H + theta[6] / (1 + P^2) - theta[7] * H
  
  PMHdt
}

param.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.439, 2.037, 17.904),
  sigma = c(0.15, 0.15, NA))

modelODE <- function(tvec, state, parameters) {
  list(as.vector(hes1modelODE(parameters, t(state), tvec)))
}

x <- deSolve::ode(y = param.true$x0, times = seq(0, 60 * 4, by = 0.01),
                  func = modelODE, parms = param.true$theta)


set.seed(12321)
y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5), ])
names(y) <- c("time", "P", "M", "H")
y$P <- y$P * exp(rnorm(nrow(y), sd = param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd = param.true$sigma[2]))


y$H <- NaN
y$P[y$time %in% seq(7.5, 240, by = 15)] <- NaN
y$M[y$time %in% seq(0, 240, by = 15)] <- NaN

pdf(file = "figures/hes1setup.pdf", width = 6, height = 4)
compnames <- c("P", "M", "H")
matplot(x[, "time"], x[, -1], type = "l", lty = 1,
        xlab = "Time (min)", ylab = "Level")
matplot(y$time, y[,-1], type = "p", col = 1:(ncol(y)-1), pch = 20, add = TRUE)
legend("topright", compnames, lty = 1, col = c("black", "red", "green"))
dev.off()

y.tilde <- y
y.tilde[, names(y.tilde) != "time"] <- log(y.tilde[, names(y.tilde) != "time"])

hes1logmodelODE <- function (theta, x, tvec) {
	P = exp(x[, 1])
	M = exp(x[, 2])
	H = exp(x[, 3])
	
	PMHdt <- array(0, c(nrow(x), ncol(x)))
	PMHdt[, 1] = -theta[1] * H + theta[2] * M / P - theta[3]
	PMHdt[, 2] = -theta[4] + theta[5] / (1 + P^2) / M
	PMHdt[, 3] = -theta[1] * P + theta[6] / (1 + P^2) / H - theta[7]

	PMHdt
}


hes1logmodelDx <- function (theta, x, tvec) {
  logP = x[, 1]
  logM = x[, 2]
  logH = x[, 3]
  
  Dx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  dP = -(1 + exp(2 * logP))^(-2) * exp(2 * logP) * 2
  Dx[, 1, 1] = -theta[2] * exp(logM - logP)
  Dx[, 2, 1] = theta[2] * exp(logM - logP)
  Dx[, 3, 1] = -theta[1] * exp(logH)
  Dx[, 1, 2] = theta[5] * exp(-logM) * dP
  Dx[, 2, 2] = -theta[5] * exp(-logM) / (1 + exp(2 * logP))
  Dx[, 1, 3] = -theta[1] * exp(logP) + theta[6] * exp(-logH) * dP
  Dx[, 3, 3] = -theta[6] * exp(-logH) / (1 + exp(2 * logP))
  	
  Dx
}


hes1logmodelDtheta <- function (theta, x, tvec) {
  logP = x[, 1]
  logM = x[, 2]
  logH = x[, 3]
  
  Dtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  Dtheta[, 1, 1] = -exp(logH)
  Dtheta[, 2, 1] = exp(logM - logP)
  Dtheta[, 3, 1] = -1
  Dtheta[, 4, 2] = -1
  Dtheta[, 5, 2] = exp(-logM) / (1 + exp(2 * logP))
  Dtheta[, 1, 3] = -exp(logP)
  Dtheta[, 6, 3] = exp(-logH) / (1 + exp(2 * logP))
  Dtheta[, 7, 3] = -1
  	
  Dtheta
}


yTest <- matrix(runif(nrow(y.tilde) * (ncol(y.tilde) - 1)),
                nrow = nrow(y.tilde), ncol = ncol(y.tilde) - 1)
thetaTest <- runif(7)
testDynamicalModel(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, 
                   "Hes1 log", yTest, thetaTest, y.tilde[, "time"])

hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = hes1logmodelDx,
  fOdeDtheta = hes1logmodelDtheta,
  thetaLowerBound = rep(0, 7),
  thetaUpperBound = rep(Inf, 7)
)

hes1result <- MagiSolver(y.tilde, hes1logmodel, control = list(sigma = param.true$sigma, useFixedSigma = TRUE))

pdf(file = "figures/hes1-trace.pdf", width = 8, height = 4)
theta.names <- c("a", "b", "c", "d", "e", "f", "g")
plot(hes1result, type = "trace", par.names = theta.names, nplotcol = 4)
dev.off()

summary(hes1result, par.names = theta.names)

pdf(file = "figures/hes1-plot-magioutput.pdf", width = 8, height = 4)
plot(hes1result, lwd = 2, col = "forestgreen", comp.names = compnames,
     xlab = "Time", ylab = "log(Level)")
dev.off()

xLB <- exp(apply(hes1result$xsampled, c(2,3),
                 function(x) quantile(x, 0.025)))
xMean <- exp(apply(hes1result$xsampled, c(2,3), mean))
xUB <- exp(apply(hes1result$xsampled, c(2,3),
                 function(x) quantile(x, 0.975)))

pdf(file = "figures/hes1-inferred-trajectories.pdf", width = 8, height = 4)
layout(rbind(c(1, 2, 3), c(4, 4, 4)), heights = c(5, 0.5))
ylim_lower <- c(1.5, 0.5, 0)
ylim_upper <- c(10.0, 3.5, 21)
compobs <- c("17 observations", "16 observations", "unobserved")
times <- y[, "time"]

for (i in 1:3) {
  plot(times, xMean[, i], type = "n", xlab = "time", ylab = compnames[i],
        ylim = c(ylim_lower[i], ylim_upper[i]))
  mtext(paste0(compnames[i], " (", compobs[i], ")"), cex = 1)
  
  polygon(c(times, rev(times)), c(xUB[, i], rev(xLB[, i])),
          col = "skyblue", border = NA)
  
  lines(x[, 1], x[, 1+i], col = "red", lwd = 2)
  lines(times, xMean[, i], col = "forestgreen", lwd = 2)
  points(times, exp(y[, 1+i]))
}

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n',
     xlab = NA, ylab = NA, frame.plot = FALSE)

legend("center", c("truth", "inferred trajectory",
                   "95% credible interval", "noisy observations"),
  lty = c(1, 1, 0, 0), lwd = c(2, 2, 0, 1), bty = "n",
  col = c("red", "forestgreen", NA, "black"), fill = c(0, 0, "skyblue", 0),
  border = c(0, 0, "skyblue", 0), pch = c(NA, NA, 15, 1), horiz = TRUE)
dev.off()


# Figure 2 (Example visualization of manifold constraint)
set.seed(1234)
dat <- data.frame(time = c(1, 3, 4, 5), y = c(-2, -1, 1, 1))

tvec <- seq(0, 10, by = 0.05)
phi <- matrix(c(10, 3))
sigma <- 0.5

## Conditional mean and covariance given observations, sample 5 trajectories
condMean <- gpmean(dat$y, dat$time, tvec, phi, sigma)
condSigma <- gpcov(dat$y, dat$time, tvec, phi, sigma)
ysim <- MASS::mvrnorm(5, condMean, condSigma)

## Set up a simple linear ODE system and use MAGI to condition on it
linearODE <- function(theta, x, tvec) {
  lindt = array(0, c(nrow(x), ncol(x)))
  lindt[, 1] = theta[1] * x[,1] + theta[2]
  lindt
}
linearDx <- function(theta, x, tvec) {
  Dx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  Dx[, 1, 1] = theta[1]
  Dx
}
linearDtheta <- function(theta, x, tvec) {
  Dtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  Dtheta[, 1, 1] = x[,1]
  Dtheta[, 2, 1] = 1
  Dtheta
}

linearmodel <- list(
  fOde = linearODE,
  fOdeDx = linearDx,
  fOdeDtheta = linearDtheta,
  thetaLowerBound = rep(-Inf, 2),
  thetaUpperBound = rep(Inf, 2)
)

yIn <- data.frame(time=tvec, y = NA)
yIn[yIn$time %in% dat$time, "y"] <- dat$y
res <- MagiSolver(yIn, linearmodel, control=list(sigma = sigma, phi = phi, useFixedSigma = TRUE, niterHmc = 5000))

## Panel (a): GP given observations only
pdf("figures/GP-given-obs.pdf", height = 4.5, width = 5)
par(mar = c(4, 4, 1.5, 1))
plot(tvec, ysim[1,], type = 'n', ylim = c(-6.5, 7.5), xlab = "Time", ylab = "x")
polygon(c(rev(tvec), tvec), c(rev(condMean + 1.96 * sqrt(diag(condSigma))), condMean - 1.96 * sqrt(diag(condSigma))), col = 'grey80', border = NA)
for (j in 1:5)
  lines(tvec, ysim[j, ], col = rainbow(12)[j+7])
points(dat$time, dat$y, pch = 20)
dev.off()

## Panel (b): GP given observations and manifold constraint
pdf("figures/GP-given-obs-constraint.pdf", height = 4.5, width = 5)
plot(res, est = "none", ci.col = "grey80", comp.names = "", ylim = c(-6.5, 7.5), xlab = "Time", ylab = "x", obs = FALSE)
for (j in 1:5)
  lines(tvec, res$xsampled[sample(1:length(res$lp), 1), , ], col = rainbow(12)[j+7])
points(dat$time, dat$y, pch = 20)
dev.off()


# Fitzhugh-Nagumo equations

data("FNdat")
set.seed(12321)

y_I0 <- setDiscretization(FNdat, by = 0.5)

y_I1 <- setDiscretization(y_I0, level = 1)
y_I2 <- setDiscretization(y_I0, level = 2)
y_I3 <- setDiscretization(y_I0, level = 3)

fnmodelODE <- function(theta, x, tvec) {
  V <- x[, 1]
  R <- x[, 2]
  
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] = theta[3] * (V - V^3 / 3.0 + R)
  result[, 2] = -1.0/theta[3] * (V - theta[1] + theta[2] * R)
  
  result
}

fnmodelDx <- function(theta, x, tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  V = x[, 1]
  
  resultDx[, 1, 1] = theta[3] * (1 - V^2)
  resultDx[, 2, 1] = theta[3]
  
  resultDx[, 1, 2] = -1.0 / theta[3]
  resultDx[, 2, 2] = -theta[2] / theta[3]
  
  resultDx
}

fnmodelDtheta <- function(theta, x, tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  V = x[, 1]
  R = x[, 2]
  
  resultDtheta[, 3, 1] = V - V^3 / 3.0 + R
  
  resultDtheta[, 1, 2] =  1.0 / theta[3]
  resultDtheta[, 2, 2] = -R / theta[3]
  resultDtheta[, 3, 2] = 1.0 / (theta[3]^2) * (V - theta[1] + theta[2] * R)
  
  resultDtheta
}

fnmodel <- list(
  fOde = fnmodelODE,
  fOdeDx = fnmodelDx,
  fOdeDtheta = fnmodelDtheta,
  thetaLowerBound = c(0, 0, 0),
  thetaUpperBound = c(Inf, Inf, Inf)
)


FNres0 <- MagiSolver(y_I0, fnmodel, control = list(niterHmc = 10000))
FNres1 <- MagiSolver(y_I1, fnmodel, control = list(niterHmc = 10000))
FNres2 <- MagiSolver(y_I2, fnmodel, control = list(niterHmc = 10000))
FNres3 <- MagiSolver(y_I3, fnmodel, control = list(niterHmc = 10000, nstepsHmc = 1000))

FNpar.names <- c("a", "b", "c", "sigmaV", "sigmaR")
FNsummary <- lapply(list(FNres0, FNres1, FNres2, FNres3),
               function(x) summary(x, sigma = TRUE, par.names = FNpar.names))

pdf(file = "figures/FNparam.pdf", width = 10, height = 3)
layout(rbind(c(1:5), rep(6, 5)), heights = c(5, 0.25))
for (i in 1:length(FNpar.names)) {
  par(mar = c(2, 4, 1.5, 1))
  estCI <- sapply(FNsummary, function(x) x[,i])
  plot(1:4, xlim = c(0, 5), ylim = c(min(estCI[2, ]), max(estCI[3, ])),
       xaxt = 'n', xlab = '', ylab = '', type = 'n')
  segments(1:4, y0 = estCI[2, ], y1 = estCI[3, ], col = 1:4, lwd = 2)
  mtext(FNpar.names[i])
  points(1:4, estCI[1, ], col = 1:4, cex = 2)
}

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n',
     xlab = NA, ylab = NA, frame.plot = FALSE)
legend("center", c("I0", "I1", "I2", "I3"),
       col = 1:4, lwd = 4, horiz = TRUE, bty = "n")
dev.off()

fnmodelODEsolve <- function(tvec, state, parameters) {
  list(as.vector(fnmodelODE(parameters, t(state), tvec)))
}

tvec <- seq(0, 20, by = 0.01)
FNcalcTraj <- function(res) {
  x0.est <- apply(res$xsampled[, 1, ], 2, mean)
  theta.est <- apply(res$theta, 2, mean)
  
  x <- deSolve::ode(y = x0.est, times = tvec,
          func = fnmodelODEsolve, parms = theta.est)
  x
}

FNtr <- lapply(list(FNres0, FNres1, FNres2, FNres3), FNcalcTraj)

pdf(file = "figures/FNtraj.pdf", width = 10, height = 5)
layout(rbind(c(1, 2), c(3, 3)), heights = c(5, 0.25))
plot(FNdat$time, FNdat$V, xlab = "Time", ylab = "V")
matplot(tvec, sapply(FNtr, function(x) x[, 2]), type = "l", lty = 1, add = TRUE)
plot(FNdat$time, FNdat$R, xlab = "Time", ylab = "R")
matplot(tvec, sapply(FNtr, function(x) x[, 3]), type = "l", lty = 1, add = TRUE)

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n', xlab = NA, ylab = NA, frame.plot = FALSE)
legend("center", c("I0", "I1", "I2", "I3"), col = 1:4, lwd = 4, horiz = TRUE, bty = "n")
dev.off()

FN.rmsd <- sapply(FNtr, function(x)
  sqrt(colMeans((subset(x, time %in% FNdat$time) - FNdat[, 2:3])^2)))

colnames(FN.rmsd) <- c("I0", "I1", "I2", "I3")
round(FN.rmsd, 3)


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

phiEst <- matrix(0, nrow = 2, ncol = ncol(y) - 1)
sigmaInit <- rep(0, ncol(y) - 1)
for (j in 1:(ncol(y) - 1)) {
  hyperparam <- gpsmoothing(y[, j+1], y[, "time"])
  phiEst[, j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
}

compnames <- c("TU", "TI", "V")
complabels <- c("Concentration", "Concentration", "Load")

pdf(file = "figures/HIVsim-fit.pdf", width = 8, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 1.5, 1))
tOut <- seq(0, 20, by = 0.025)
for (j in 1:3) {
  plot(y[, "time"], y[, j+1], type = 'n',
      xlab = "Time", ylab = complabels[j])
  mtext(compnames[j])
  fitMean <- gpmean(y[, j+1], y[, "time"], tOut, phiEst[, j], sigmaInit[j])
  fitCov <- gpcov(y[, j+1], y[, "time"], tOut, phiEst[, j], sigmaInit[j])
  gp_UB <- fitMean + 1.96 * sqrt(diag(fitCov))
  gp_LB <- fitMean - 1.96 * sqrt(diag(fitCov))

  polygon(c(tOut, rev(tOut)), c(gp_UB, rev(gp_LB)),
          col = "grey80", border = NA)
  points(y[, "time"], y[, j+1], cex = 0.5)
  lines(tOut, fitMean, col = "steelblue")
}

colnames(phiEst) <- compnames
phiEst
sigmaInit

phiEst[, 3] <- c(1e7, 0.5)
sigmaInit[3] <- 100

j <- 3
fitMean <- gpmean(y[, j+1], y[, "time"], tOut, phiEst[, j], sigmaInit[j])
fitCov <- gpcov(y[, j+1], y[, "time"], tOut, phiEst[, j], sigmaInit[j])
gp_UB <- fitMean + 1.96 * sqrt(diag(fitCov))
gp_LB <- fitMean - 1.96 * sqrt(diag(fitCov))

polygon(c(tOut, rev(tOut)), c(gp_UB, rev(gp_LB)),
        col = "grey80", border = NA)
lines(tOut, fitMean, col = "darkorange")

dev.off()

y_I <- setDiscretization(y, level = 1)
HIVresult <- MagiSolver(y_I, hivtdmodel,
                  control = list(phi = phiEst, sigma = sigmaInit))

summary(HIVresult, sigma = TRUE, par.names =
   c("lambda", "rho", "delta", "N", "c", "sigma_TU", "sigma_TI", "sigma_V"))

xMean <- apply(HIVresult$xsampled, c(2, 3), mean)
xLB <- apply(HIVresult$xsampled, c(2, 3), function(x) quantile(x, 0.025))
xUB <- apply(HIVresult$xsampled, c(2, 3), function(x) quantile(x, 0.975))

pdf(file = "figures/HIVinf.pdf", width = 8, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 1.5, 1))
for (i in 1:3) {
  plot(y_I$time, xMean[, i], type = "n", xlab = "Time", ylab = complabels[i])
  mtext(compnames[i])
  polygon(c(y_I$time, rev(y_I$time)), c(xUB[, i], rev(xLB[, i])),
          col = "skyblue", border = NA)  
  lines(y_I$time, xMean[, i], col = "forestgreen", lwd = 2)
  lines(param.true$times, xtrue[, i+1], col = "red", lwd = 1)
}
dev.off()


# Hamiltonian Monte Carlo - example of sticky samples with high autocorrelation

data("FNdat")
set.seed(12321)

y_I0 <- setDiscretization(FNdat, by = 0.5)
y_I3 <- setDiscretization(y_I0, level = 3)
FNres3b <- MagiSolver(y_I3, fnmodel, control = list(niterHmc = 10000))

pdf(file = "figures/fn-sticky-trace.pdf", width = 8, height = 4)
plot(FNres3b, type = "trace", par.names = FNpar.names, sigma = TRUE)
dev.off()

sessionInfo()
