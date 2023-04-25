#' get posterior mean curve for value conditioning on observed y, phi, sigma
#'
#' use to visulize Gaussian process smoothing without ODE, and for initialize
#' the x theta sampler
#'
#' @noRd
getMeanCurve <- function(x, y, x.new, phi.mat, sigma.mat, kerneltype="matern", deriv=FALSE){
  tvec <- c(x.new,x)

  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2

  signr <- -sign(foo)

  y.new <- matrix(NA, nrow(phi.mat), length(x.new))
  dy.new <- matrix(NA, nrow(phi.mat), length(x.new))

  for(it in 1:nrow(phi.mat)){
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]

    covObj <- calCov(phi, r, signr, complexity = as.numeric(deriv), kerneltype=kerneltype)
    C <- covObj$C

    diag(C)[-(1:length(x.new))] <- diag(C)[-(1:length(x.new))]+sigma^2
    y.new[it, ] <- C[1:length(x.new),-(1:length(x.new))]%*%solve(C[-(1:length(x.new)),-(1:length(x.new))], y)
    if(deriv){
      dy.new[it, ] <- covObj$Cprime[1:length(x.new), 1:length(x.new)] %*% solve(covObj$C[1:length(x.new), 1:length(x.new)], y.new[it, ])
    }
  }

  if(deriv){
    return(list(y.new, dy.new))
  }else{
    return(y.new)
  }
}

#' Conditional mean of Gaussian process given observations
#' 
#' Compute the conditional mean of a Gaussian process (and optionally, its derivative), given a vector of observations, hyper-parameters \code{phi}, and noise standard deviation \code{sigma}.
#' 
#' @param yobs vector of observations
#' @param tvec vector of time points corresponding to observations
#' @param tnew vector of time points at which the conditional mean should be computed
#' @param phi vector of hyper-parameters for the covariance kernel (\code{kerneltype})
#' @param sigma noise standard deviation of the observations
#' @param kerneltype the covariance kernel, types \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}, \code{generalMatern} are supported.  See \code{\link{calCov}} for their definitions.
#' @param sigma the noise level (if known). By default, both \code{phi} and \code{sigma} are estimated. If a value for \code{sigma} is supplied, then \code{sigma} is held fixed at the supplied value and only \code{phi} is estimated.
#' @param deriv logical; if true, the conditional mean of the GP's derivative is also computed
#' 
#' @return A vector with the values of the conditional mean function evaluated at the time points in \code{tnew}. If \code{deriv = TRUE}, returned with an additional attribute \code{deriv} that contains the values of the conditional mean of the GP derivative evaluated at the time points in \code{tnew}.
#'
#' @examples 
#' # Load Fitzhugh-Nagumo dataset
#' data(FNdat)
#' 
#' tnew <- seq(0, 20, by = 0.5)
#' 
#' # GP mean of V component at time points in tnew given observations
#' gpmean(FNdat$V, FNdat$time, tnew, c(2.3, 1.2), 0.2)
#' 
#' @export
gpmean <- function(yobs, tvec, tnew, phi, sigma, kerneltype="generalMatern", deriv=FALSE) {
  res <- getMeanCurve(tvec, yobs, tnew, as.matrix(t(phi)), sigma, kerneltype, deriv)
  
  if (deriv) {
    ret <- res[[1]]
    attr(ret, "deriv") <- res[[2]]
  } else {
    ret <- res
  }
  
  ret
}

#' Conditional covariance of Gaussian process given observations
#' 
#' Compute the conditional covariance of a Gaussian process, given a vector of observations, hyper-parameters \code{phi}, and noise standard deviation \code{sigma}.
#' 
#' @param yobs vector of observations
#' @param tvec vector of time points corresponding to observations
#' @param tnew vector of time points at which the conditional covariance should be computed
#' @param phi vector of hyper-parameters for the covariance kernel (\code{kerneltype})
#' @param sigma noise standard deviation of the observations
#' @param kerneltype the covariance kernel, types \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}, \code{generalMatern} are supported.  See \code{\link{calCov}} for their definitions.
#' @param sigma the noise level (if known). By default, both \code{phi} and \code{sigma} are estimated. If a value for \code{sigma} is supplied, then \code{sigma} is held fixed at the supplied value and only \code{phi} is estimated.
#' 
#' @return The conditional covariance matrix for the GP evaluated at the time points in \code{tnew}.
#'
#' @examples 
#' # Load Fitzhugh-Nagumo dataset
#' data(FNdat)
#' 
#' tnew <- seq(15, 20, by = 0.5)
#' 
#' # GP covariance of V component at time points in tnew given observations
#' gpcov(FNdat$V, FNdat$time, tnew, c(2.3, 1.2), 0.2)
#' 
#' @export
gpcov <- function(yobs, tvec, tnew, phi, sigma, kerneltype="generalMatern") {

  tvec <- c(tnew, tvec)
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  signr <- -sign(foo)

  covObj <- calCov(phi, r, signr, complexity = 0, kerneltype=kerneltype)
  C <- covObj$C
    
  diag(C)[-(1:length(tnew))] <- diag(C)[-(1:length(tnew))]+sigma^2
  ret <- C[1:length(tnew),1:length(tnew)] - C[1:length(tnew),-(1:length(tnew))] %*% solve(C[-(1:length(tnew)),-(1:length(tnew))]) %*% t(C[1:length(tnew),-(1:length(tnew))])
  0.5 * (ret + t(ret)) # ensure symmetric
}



#' Set discretization level
#'
#' @description Set the discretization level of a data matrix for input to \code{\link{MagiSolver}}, by inserting time points where the GP is constrained to the derivatives of the ODE system.
#' 
#' @param dat data matrix. Must include a column with name `time`.
#' @param level discretization level (a positive integer). \code{2^level - 1} equally-spaced time points will be inserted between each row of \code{dat}.
#' @param by discretization interval. As an alternative to \code{level}, time points will be inserted (as needed) to form an equally-spaced discretization set from the first to last observations of \code{dat}, with interval \code{by} between successive discretization points. This can be useful when the time points in \code{dat} are unevenly spaced.
#'
#' @details 
#' Specify the desired discretization using \code{level} or \code{by}.
#' 
#' @return Returns a data matrix with the same columns as \code{dat}, with rows added for the inserted discretization time points.
#' 
#' @examples
#' dat <- data.frame(time = 0:10, x = rnorm(11))
#' setDiscretization(dat, level = 2)
#' setDiscretization(dat, by = 0.2)
#'
#' @export
setDiscretization <- function(dat, level, by) {
  if (!("time" %in% colnames(dat)))
    stop("\"dat\" is missing a column \"time\"")
  if (ncol(dat) < 2)
    stop("\"dat\" does not have any components")

  if (!missing(level) & !missing(by))
    stop("specify either \"level\" or \"by\" but not both")
  if (missing(level) & missing(by)) {
    warning("no discretization added")
    level = 0
  }
  
  if (!missing(level)) {
    if (round(level) >= 0) {
      if(round(level) == 0){
        return(dat)
      }
      newdata <- dat
      newdata <- newdata[order(newdata$time),]
      dumdat <- newdata[-1,]
      dumdat[] <- NaN
      dumdat$time <- (newdata$time[-1] + newdata$time[-nrow(newdata)])/2
      newdata <- rbind(newdata, dumdat)
      newdata <- newdata[order(newdata$time),]
      return(setDiscretization(newdata, level-1))
  
    } else {
      stop("\"level\" must be a non-negative integer")
    }
  }
  
  if (!missing(by)) {
    fillC <- seq(min(dat[,"time"]), max(dat[,"time"]), by = by)
    newdata <- data.frame(time = sort(unique(c(fillC, dat[,"time"]))))
    newdata <- cbind(newdata, matrix(NaN, nrow = length(newdata$time), ncol = ncol(dat)-1 ))
    datd <- dat[,!colnames(dat)=="time",drop=FALSE]
    for (i in 1:nrow(newdata)) {
      loc <- match( newdata$time[i], dat[, "time"])
      if (!is.na(loc))
        newdata[i,2:ncol(dat)] = datd[loc,]
    }
    colnames(newdata) <- c("time",colnames(datd))
    
    return(newdata)
  }
}
