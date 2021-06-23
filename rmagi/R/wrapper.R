hes1modelODE <- function(theta, x, tvec) {
  P = x[,1]
  M = x[,2]
  H = x[,3]

  PMHdt = array(0, c(nrow(x), ncol(x)))
  PMHdt[,1] = -theta[1]*P*H + theta[2]*M - theta[3]*P
  PMHdt[,2] = -theta[4]*M + theta[5]/(1+P^2)
  PMHdt[,3] = -theta[1]*P*H + theta[6]/(1+P^2) - theta[7]*H

  PMHdt
}

param.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.438575, 2.037488, 17.90385),
  sigma = c(0.15, 0.15, NA)
)

modelODE <- function(t, state, parameters) {
  list(as.vector(hes1modelODE(parameters, t(state), t)))
}
x <- deSolve::ode(y = param.true$x0, times = seq(0, 60*4, by = 0.01), func = modelODE, parms = param.true$theta)

y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5),])
names(y) <- c("time", "P", "M", "H")

set.seed(1234)
y$P <- y$P * exp(rnorm(nrow(y), sd=param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd=param.true$sigma[2]))

y$H <- NaN
y$P[y$time %in% seq(7.5,240,by=15)] <- NaN
y$M[y$time %in% seq(0,240,by=15)] <- NaN

pdf(file="../Fig1.pdf", width=6, height=4)
matplot(x[, "time"], x[, -1], type="l", lty=1, xlab="Time (min)", ylab="Level")
matplot(y$time, y[,-1], type="p", col=1:(ncol(y)-1), pch=20, add = TRUE)
legend("topright", c("P", "M", "H"), lty=1, col=c("black", "red", "green"))
dev.off()

library(magi)


### Temp (simple) R wrapper
MagiSolver <- function(y, odeModel, tvec, control = list()) {

  if (!is.null(control$sigma))
    sigmaExogenous = control$sigma
  else
    sigmaExogenous = numeric(0)

  if (!is.null(control$phi))
    phiExogenous = control$phi
  else
    phiExogenous = matrix(numeric(0))

  if (!is.null(control$xInit))
    xInitExogenous = control$xInit
  else
    xInitExogenous = matrix(numeric(0))

  if (!is.null(control$thetaInit))
    thetaInitExogenous = control$thetaInit
  else
    thetaInitExogenous = matrix(numeric(0))

  if (!is.null(control$mu))
    muExogenous = control$mu
  else
    muExogenous = matrix(numeric(0))

  if (!is.null(control$dotmu))
    dotmuExogenous = control$dotmu
  else
    dotmuExogenous = matrix(numeric(0))

  if (!is.null(control$priorTemperature)) {
    priorTemperatureLevel = control$priorTemperature
    priorTemperatureDeriv = control$priorTemperature
  } else {
    priorTemperatureLevel = 1/mean(!is.na(y))
    priorTemperatureDeriv = 1/mean(!is.na(y))
  }

  if (!is.null(control$niterHmc))
    niterHmc = control$niterHmc
  else
    niterHmc = 20000

  if (!is.null(control$burninRatio))
    burninRatio = control$burninRatio
  else
    burninRatio = 0.5

  if (!is.null(control$nstepsHmc))
    nstepsHmc = control$nstepsHmc
  else
    nstepsHmc = 100

  if (!is.null(control$stepSizeFactor))
    stepSizeFactor = control$stepSizeFactor
  else
    stepSizeFactor = 0.01

  if (!is.null(control$bandSize))
    bandSize = control$bandSize
  else
    bandSize = 20

  if (!is.null(control$useFixedSigma))
    useFixedSigma = control$useFixedSigma
  else
    useFixedSigma = FALSE


  samplesCpp <- magi:::solveMagiRcpp(
    yFull = data.matrix(y),
    odeModel = odeModel,
    tvecFull = tvec,
    sigmaExogenous = sigmaExogenous,
    phiExogenous = phiExogenous,
    xInitExogenous = xInitExogenous,
    thetaInitExogenous = thetaInitExogenous,
    muExogenous = muExogenous,
    dotmuExogenous = dotmuExogenous,
    priorTemperatureLevel = priorTemperatureLevel,
    priorTemperatureDeriv = priorTemperatureDeriv,
    priorTemperatureObs = 1,
    kernel = "generalMatern",
    nstepsHmc = nstepsHmc,
    burninRatioHmc = burninRatio,
    niterHmc = niterHmc,
    stepSizeFactorHmc = stepSizeFactor,
    nEpoch = 1,
    bandSize = bandSize,
    useFrequencyBasedPrior = TRUE,
    useBand = TRUE,
    useMean = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = useFixedSigma,
    verbose = TRUE)

  phiUsed <- samplesCpp$phi
  samplesCpp <- samplesCpp$llikxthetasigmaSamples

  samplesCpp <- samplesCpp[,,1]

  out <- samplesCpp[-1,1,drop=FALSE]

  llikId <- 1
  xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(y)))
  thetaId <- (max(xId)+1):(max(xId)+length(odeModel$thetaLowerBound))
  sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(y))

  burnin <- as.integer(niterHmc*burninRatio)
  return(
    list(theta=t(samplesCpp[thetaId, -(1:burnin)]),
         xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                        dim=c(niterHmc-burnin, nrow(y), ncol(y))),
         lp=samplesCpp[llikId,-(1:burnin)],
         sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]),
         phi = phiUsed)
  )


}

hes1logmodelODE <- function (theta, x, tvec) {
  eP = exp(x[, 1])
  eM = exp(x[, 2])
  eH = exp(x[, 3])

  PMHdt <- array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * eH + theta[2] * eM/eP - theta[3]
  PMHdt[, 2] = -theta[4] + theta[5]/(1 + eP^2)/eM
  PMHdt[, 3] = -theta[1] * eP + theta[6]/(1 + eP^2)/eH - theta[7]
  PMHdt
}

hes1logmodelDx <- function (theta, x, tvec) {
  P = x[, 1]
  M = x[, 2]
  H = x[, 3]

  Dx <- array(0, c(nrow(x), ncol(x), ncol(x)))

  dP = -(1 + exp(2 * P))^(-2) * exp(2 * P) * 2
  Dx[, 1, 1] = -theta[2] * exp(M - P)
  Dx[, 2, 1] = theta[2] * exp(M - P)
  Dx[, 3, 1] = -theta[1] * exp(H)
  Dx[, 1, 2] = theta[5] * exp(-M) * dP
  Dx[, 2, 2] = -theta[5] * exp(-M)/(1 + exp(2 * P))
  Dx[, 1, 3] = -theta[1] * exp(P) + theta[6] * exp(-H) * dP
  Dx[, 3, 3] = -theta[6] * exp(-H)/(1 + exp(2 * P))

  Dx
}

hes1logmodelDtheta <- function (theta, x, tvec) {

  P = x[, 1]
  M = x[, 2]
  H = x[, 3]

  Dtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  Dtheta[, 1, 1] = -exp(H)
  Dtheta[, 2, 1] = exp(M - P)
  Dtheta[, 3, 1] = -1
  Dtheta[, 4, 2] = -1
  Dtheta[, 5, 2] = exp(-M)/(1 + exp(2 * P))
  Dtheta[, 1, 3] = -exp(P)
  Dtheta[, 6, 3] = exp(-H)/(1 + exp(2 * P))
  Dtheta[, 7, 3] = -1

  Dtheta
}

hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = hes1logmodelDx,
  fOdeDtheta = hes1logmodelDtheta,
  thetaLowerBound = rep(0,7),
  thetaUpperBound = rep(Inf,7)
)

result <- MagiSolver(log(y[,-1]), hes1logmodel, y$time, control=list(sigma = c(0.15,0.15,0.00001), useFixedSigma = TRUE))

