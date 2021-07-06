#' MAnifold-constrained Gaussian process Inference (MAGI)
#'
#' @description Main function for the MAGI method for inferring the parameters and trajectories of dynamic systems governed by ordinary differential equations.
#' 
#' @param y data matrix of observations
#' @param odeModel list of ODE functions
#' @param tvec vector of discretization time points
#' @param control list of control variables, including `sigma`, `phi`, `xInit`, `thetaInit`, `mu`, `dotmu`, `priorTemperature`, `niterHmc`
#' `burninRatio`, `nstepsHmc`, `stepSizeFactor`, `bandSize`, `useFixedSigma`.  See details.
#' 
#' @details 
#' \describe{
#'   \item{\code{sigma}}{First item}
#'   \item{\code{phi}}{Second item}
#' }
#' 
#' @export 
MagiSolver <- function(y, odeModel, tvec, control = list()) {
  if (missing(tvec) & is.null(y$time))
    stop("must supply a \"time\" column in y or specify a vector of time points via the argument \"tvec\"")
  
  if (missing(tvec)) {
    tvec <- y$time
    y <- y[, !names(y)=="time", drop=FALSE]
  }
  
  if (!is.null(control$sigma)) {
    sigmaExogenous = control$sigma
    sigmaExogenous[!is.finite(sigmaExogenous)] <- 1
  }
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
    nstepsHmc = 200

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


  samplesCpp <- solveMagiRcpp(
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
