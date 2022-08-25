# testhat helper functions

#' full log likelihood in R with contribution from each component
#' 
#' value is exact for x, theta, phi, sigma. we simply need to sample from this
#' if computation is feasible. phi info is supplied explicitly
#' 
#' @noRd
loglikWithNormalizingConstants <- function(x, theta, phi, sigma, y, rInput, signrInput, kerneltype = "matern")  {
  CovV <- calCov(phi[1:2], rInput, signrInput, complexity = 3, kerneltype=kerneltype)
  CovR <- calCov(phi[3:4], rInput, signrInput, complexity = 3, kerneltype=kerneltype)
  lambda <- 1
  
  if(length(lambda) < 3) lambda <- c(lambda, rep(1, 3-length(lambda)))
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) { return(1e9)}
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- sum(is.finite(y[,1]))
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( sum(dnorm(y[is.finite(y[,1]),1], Vsm[is.finite(y[,1])], sigma, log=TRUE)) * lambda[1], 
                mvtnorm::dmvnorm(t(fr), sigma = CovV$Kphi, log=TRUE) * lambda[2],  
                mvtnorm::dmvnorm(t(Vsm), sigma = CovV$C, log=TRUE) * lambda[3])
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( sum(dnorm(y[is.finite(y[,2]),2], Rsm[is.finite(y[,2])], sigma, log=TRUE)) * lambda[1], 
                mvtnorm::dmvnorm(t(fr), sigma = CovR$Kphi, log=TRUE) * lambda[2],  
                mvtnorm::dmvnorm(t(Rsm), sigma = CovR$C, log=TRUE) * lambda[3])
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
}

