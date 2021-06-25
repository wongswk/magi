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

#' insert nan in simulated data for explicit control of discretization
#'
#' @param mydata a data frame that contains at least one column `time`
#' @param level the level to insert nan: 2^level - 1 points will be inserted between two data points
#'
#' @export
insertNaN <- function(mydata, level){
  if(level==0){
    return(mydata)
  }
  newdata <- mydata
  newdata <- newdata[order(newdata$time),]
  dummydata <- newdata[-1,]
  dummydata[] <- NaN
  dummydata$time <- (newdata$time[-1] + newdata$time[-nrow(newdata)])/2
  newdata <- rbind(newdata, dummydata)
  newdata <- newdata[order(newdata$time),]
  return(insertNaN(newdata, level-1))
}
