#' `magi`: MAnifold-Constrained Gaussian Process Inference
#'
#' `magi` is a package that provides fast and accurate inference for the parameter estimation problem in Ordinary Differential Equations, including the case when there are unobserved system components.
#' In the references below, please see Yang, Wong, and Kou (2021) for details of the MAGI method (MAnifold-constrained Gaussian process Inference), and Wong, Yang, and Kou (2022) for a detailed user guide.
#'
#' @references
#'
#' Yang, S., Wong, S. W. K., & Kou, S. C. (2021). Inference of Dynamic Systems from Noisy and 
#' Sparse Data via Manifold-constrained Gaussian Processes. *Proceedings of the National Academy of Sciences*, 118 (15), e2020397118. \doi{10.1073/pnas.2020397118}
#'
#' Wong, S. W. K., Yang, S., & Kou, S. C. (2022). `MAGI`: A Package for Inference of Dynamic Systems from Noisy and Sparse Data via Manifold-constrained Gaussian Processes. \url{https://arxiv.org/abs/2203.06066}
#'
#' @name magi
#' @md
NULL
#> NULL


#' MagiSolver output (\code{magioutput}) object 
#' @description Check for and create a magioutput object
#' @param object an R object
#' @param ... arguments required to create a magioutput object. See details.
#' @return logical. Is the input a magioutput object?
#' @details
#' Using the core \code{\link{MagiSolver}} function returns a \code{magioutput} object as output, which is a list that contains the following elements:
#' \describe{
#' \item{\code{theta}}{matrix of MCMC samples for the system parameters \eqn{\theta}, after burn-in.}
#' \item{\code{xsampled}}{array of MCMC samples for the system trajectories at each discretization time point, after burn-in.}
#' \item{\code{sigma}}{matrix of MCMC samples for the observation noise SDs \eqn{\sigma}, after burn-in.}
#' \item{\code{phi}}{matrix of estimated GP hyper-parameters, one column for each system component.}
#' \item{\code{lp}}{vector of log-posterior values at each MCMC iteration, after burn-in.}
#' } 
#' Printing a \code{magioutput} object displays a brief summary of the settings used for the \code{MagiSolver} run.
#' The summary method for a \code{magioutput} object prints a table of parameter estimates, see \code{\link{summary.magioutput}} for more details. 
#' Plotting a \code{magioutput} object shows the inferred trajectories for each component, see \code{\link{plot.magioutput}} for more details.
#' 
#' 
#' @examples
#' # Set up odeModel list for the Fitzhugh-Nagumo equations
#' fnmodel <- list(
#'   fOde = fnmodelODE,
#'   fOdeDx = fnmodelDx,
#'   fOdeDtheta = fnmodelDtheta,
#'   thetaLowerBound = c(0, 0, 0),
#'   thetaUpperBound = c(Inf, Inf, Inf)
#' )
#'
#' # Example FN data
#' data(FNdat)
#'
#' # Create magioutput from a short MagiSolver run (demo only, more iterations needed for convergence)
#' result <- MagiSolver(FNdat, fnmodel, control = list(nstepsHmc = 5, niterHmc = 50)) 
#' 
#' is.magioutput(result)
#' 
#' @export

is.magioutput <- function(object) {
  inherits(object, "magioutput")
}

#' @rdname is.magioutput
#' @export
magioutput <- function(...) {
  ell <- list(...)
  class(ell) <- "magioutput"
  ell
}


#' @export
print.magioutput <- function(x, ...) {
  
  output <- paste0("MagiSolver fitted for ODE system with ", dim(x$xsampled)[3], " components and ", ncol(x$theta), " system parameters")

  output[2] <- paste0("Number of MCMC samples after burn-in = ", length(x$lp))
  output[3] <- paste0("Number of discretization points = ", dim(x$xsampled)[2])
  
  cat(output, sep = "\n")
}


#' Summary of parameter estimates from \code{magioutput} object
#' @description Computes a summary table of parameter estimates from the output of \code{MagiSolver}
#' @param object a \code{magioutput} object.
#' @param sigma logical; if true, the noise levels \eqn{\sigma} will be included in the summary.
#' @param par.names vector of parameter names for the summary table. If provided, should be the same length as the number of parameters in \eqn{\theta}, or the combined length of \eqn{\theta} and \eqn{\sigma} when \code{sigma = TRUE}.
#' @param est string specifying the posterior quantity to treat as the estimate. Default is "mean", which treats the posterior mean as the estimate. Can be one of "mean", "median", or "mode". 
#' @param lower the lower quantile of the credible interval, default is 0.025.
#' @param upper the upper quantile of the credible interval, default is 0.975.
#' @param digits integer; the number of significant digits to print.
#' @param ... additional arguments affecting the summary produced.
#' @return Returns a matrix where rows display the posterior mean, lower credible limit, and upper credible limit of each parameter.
#' @details
#' Computes parameter estimates and credible intervals from the MCMC samples. By default, the posterior mean is treated as the parameter estimate, and \code{lower = 0.025} and \code{upper = 0.975} produces a central 95\% credible interval.
#' 
#' @examples
#' # Set up odeModel list for the Fitzhugh-Nagumo equations
#' fnmodel <- list(
#'   fOde = fnmodelODE,
#'   fOdeDx = fnmodelDx,
#'   fOdeDtheta = fnmodelDtheta,
#'   thetaLowerBound = c(0, 0, 0),
#'   thetaUpperBound = c(Inf, Inf, Inf)
#' )
#'
#' # Example FN data
#' data(FNdat)
#'
#' # Create magioutput from a short MagiSolver run (demo only, more iterations needed for convergence)
#' result <- MagiSolver(FNdat, fnmodel, control = list(nstepsHmc = 5, niterHmc = 100)) 
#' 
#' summary(result, sigma = TRUE, par.names = c("a", "b", "c", "sigmaV", "sigmaR"))
#' @export
summary.magioutput <- function(object, sigma = FALSE, par.names, est = "mean", lower = 0.025, upper = 0.975, digits = 3, ...) {

  if (!is.magioutput(object)) 
    stop("\"object\" must be a magioutput object")
  
  if (est == "mean") {
    f <- mean
    est.lab <- "Mean"
  }
  if (est == "median") {
    f <- median
    est.lab <- "50%"
  }
  if (est == "mode") {
    lpmaxInd = which.max(object$lp)
    f <- function(x) x[lpmaxInd]
    est.lab <- "Mode"
  }
    
  
  theta.est <- apply(object$theta, 2,
                     function(x) c(f(x), quantile(x, lower), quantile(x, upper)))
  
  if (sigma) {
    sigma.est <- apply(object$sigma, 2,
                       function(x) c(f(x), quantile(x, lower), quantile(x, upper)))
    theta.est <- cbind(theta.est, sigma.est)
  }
  
  if (missing(par.names)) {
    par.names = paste0("theta[", 1:ncol(object$theta), "]")
    if (sigma)
      par.names = c(par.names, paste0("sigma[", 1:ncol(object$sigma), "]"))
  } else if (length(par.names) != ncol(object$theta) + sigma * ncol(object$sigma)) {
    stop(paste("vector of par.names should be length", ncol(object$theta) + sigma * ncol(object$sigma), "to match the number of parameters"))
  }
  
  colnames(theta.est) <- par.names
  rownames(theta.est) <- c(est.lab, paste0(lower*100, "%"), paste0(upper*100, "%"))
  signif(theta.est, digits)
  
}
