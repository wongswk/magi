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

#' Set discretization level
#'
#' @description Set the discretization level of a data matrix for input to \code{\link{MagiSolver}}, by inserting time points where the GP is constrained to the derivatives of the ODE system.
#' 
#' @param dat data matrix. Must include a column with name `time`.
#' @param level discretization level (a positive integer). \code{2^level - 1} equally-spaced points will be inserted between existing data points in \code{dat}.
#' @param by discretization interval. As an alternative to \code{level}, equally-spaced spaced time points will be inserted with interval \code{by} between successive points.
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
  if (is.null(dat[,"time"]))
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
