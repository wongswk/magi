#' Gaussian process smoothing
#'
#' Estimate hyper-parameters \code{phi} and noise standard deviation \code{sigma} for a vector of observations using Gaussian process smoothing.
#' 
#' @param yobs vector of observations
#' @param tvec vector of time points corresponding to observations
#' @param kerneltype the covariance kernel, types \code{matern}, \code{compact1}, \code{periodicMatern}, \code{generalMatern} are supported.  See \code{\link{calCov}} for their definitions.
#' @param sigma the noise level (if known). By default, both \code{phi} and \code{sigma} are estimated. If a value for \code{sigma} is supplied, then \code{sigma} is held fixed at the supplied value and only \code{phi} is estimated.
#' 
#' @return A list containing the elements \code{phi} and \code{sigma} with their estimated values.
#'
#' @examples
#' # Sample data and observation times
#' tvec <- seq(0, 20, by = 0.5)
#' y <- c(-1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, 
#' 0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, 
#' -0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, 
#' -1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84)
#'
#' gpsmoothing(y, tvec)
#'  
#'
#' @export
gpsmoothing <- function(yobs, tvec, kerneltype = "generalMatern", sigma = NULL) {
  
  distInput  <- abs(outer(tvec, t(tvec),'-')[,1,])
  yInput <- data.matrix(yobs - mean(yobs))
  ret <- list()
  
  if (is.null(sigma)) {
    res <- gpsmooth(yInput, distInput, kerneltype, sigmaExogenScalar = -1, TRUE)
    ret$sigma <- tail(res, 1)
    ret$phi <- res[-length(res)]
  } else {
    res <- gpsmooth(yInput, distInput, kerneltype, sigmaExogenScalar = sigma, TRUE)
    ret$sigma <- sigma
    ret$phi <- res
  }
  
  ret
}
