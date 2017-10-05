
loglikVmis <- function(x, theta, CovV, CovR, sigma, y)  {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) { return(1e9)}
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( 0, -0.5 * as.numeric(determinant(CovV$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovR$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}


loglikRmis <- function(x, theta, CovV, CovR, sigma, y)  {
  # R mis don't use 2nd col of Y
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) { return(1e9)}
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovV$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( 0, -0.5 * as.numeric(determinant(CovR$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}

xthetallikRmis <- function(x, theta, CovV, CovR, sigma, y, grad = F)  {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) {
    ret <- -1e+9
    attr(ret,"grad") <- rep(1e9, (nobs*2+3))
    return(ret)
    
  }
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  frV <- (f[,1] - CovV$mphi %*% Vsm)
  
  # R
  #CovR <- calCov(phi[3:4])
  frR <- (f[,2] - CovR$mphi %*% Rsm)
  
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2, -0.5 * t(frV) %*% CovV$Kinv %*% frV,   - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  res[2,] <- c( 0, -0.5 * t(frR) %*% CovR$Kinv %*% frR,   - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  
  if(grad) {
    # V contrib
    Vtemp <- diag( c*(1 - x[,1]^2)) - CovV$mphi
    Rtemp <- diag( rep(c,n))
    aTemp <- rep(0,n)
    bTemp <- rep(0,n)
    cTemp <- f[,1] / c
    VC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovV$Kinv %*% frV
    
    # R contrib
    Vtemp <- diag( rep( -1/c, n) )
    Rtemp <- diag( rep( -b/c, n) ) - CovR$mphi
    aTemp <- rep(1/c,n)
    bTemp <- -Rsm/c
    cTemp <- f[,2] * (-1/c)
    RC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovR$Kinv %*% frR
    
    C3 <- c(2 * CovV$Cinv %*% Vsm,  2 * CovR$Cinv %*% Rsm ,0,0,0)
    C1 <- c( 2 * (Vsm - y[,1]) / sigma^2 ,  rep(0,nobs) ,0,0,0)
    
    attr(ret,"grad") <- c((VC2 + RC2 + C3 + C1) * (-0.5))
  }
  
  return(ret)
  
  
}

xthetallikVmis <- function(x, theta, CovV, CovR, sigma, y, grad = F)  {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) {
    ret <- -1e+9
    attr(ret,"grad") <- rep(1e9, (nobs*2+3))
    return(ret)
    
  }
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  frV <- (f[,1] - CovV$mphi %*% Vsm)
  
  # R
  #CovR <- calCov(phi[3:4])
  frR <- (f[,2] - CovR$mphi %*% Rsm)
  
  res[1,] <- c( 0, -0.5 * t(frV) %*% CovV$Kinv %*% frV,   - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2, -0.5 * t(frR) %*% CovR$Kinv %*% frR,   - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  
  if(grad) {
    # V contrib
    Vtemp <- diag( c*(1 - x[,1]^2)) - CovV$mphi
    Rtemp <- diag( rep(c,n))
    aTemp <- rep(0,n)
    bTemp <- rep(0,n)
    cTemp <- f[,1] / c
    VC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovV$Kinv %*% frV
    
    # R contrib
    Vtemp <- diag( rep( -1/c, n) )
    Rtemp <- diag( rep( -b/c, n) ) - CovR$mphi
    aTemp <- rep(1/c,n)
    bTemp <- -Rsm/c
    cTemp <- f[,2] * (-1/c)
    RC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovR$Kinv %*% frR
    
    C3 <- c(2 * CovV$Cinv %*% Vsm,  2 * CovR$Cinv %*% Rsm ,0,0,0)
    C1 <- c( rep(0,nobs),  2 * (Rsm - y[,2]) / sigma^2,0,0,0)
    
    attr(ret,"grad") <- c((VC2 + RC2 + C3 + C1) * (-0.5))
  }
  
  return(ret)
  
  
}

xthU <- function(q, grad=FALSE, lambda=1) {
  x <- cbind(q[1:nobs], q[(nobs+1):(nobs*2)])
  theta <- q[(nobs*2+1):(nobs*2+3)]
  
  #xthetallik(x,theta, c(1.9840824, 1.1185157, 0.9486433, 3.268243), 0.1, fn.sim[,1:2], grad)
  xthetallik(x,theta, curCovV, curCovR, cursigma, fn.sim[,1:2], grad, lambda)
}

xthUnoODE <- function(q, grad=FALSE) {
  x <- cbind(q[1:nobs], q[(nobs+1):(nobs*2)])
  
  xthetalliknoODE(x, curCovV, curCovR, cursigma, fn.sim[,1:2], grad)
}

xthURmis <- function(q, grad=FALSE) {
  x <- cbind(q[1:nobs], q[(nobs+1):(nobs*2)])
  theta <- q[(nobs*2+1):(nobs*2+3)]
  
  xthetallikRmis(x,theta, curCovV, curCovR, cursigma, fn.sim[,1:2], grad)
}

xthUVmis <- function(q, grad=FALSE) {
  x <- cbind(q[1:nobs], q[(nobs+1):(nobs*2)])
  theta <- q[(nobs*2+1):(nobs*2+3)]
  
  xthetallikVmis(x,theta, curCovV, curCovR, cursigma, fn.sim[,1:2], grad)
}


#' marginal log likelihood without ODE on theta
#' 
#' used for Gaussian process smoothing, same as phisigllik
#' 
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
logliknoODE.mar <- function(CovV, CovR, sigma, y)  {
  n <- nrow(y)
  
  #f <- fODE(theta, x)
  res <- c(0,0)
  
  # V 
  #CovV <- calCov(phi[1:2])
  #fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1] <- dmvnorm(y[,1], sigma = CovV$C+diag(sigma^2, nrow = n), log=TRUE)
  # R
  #CovR <- calCov(phi[3:4])
  #fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2] <- dmvnorm(y[,2], sigma = CovR$C+diag(sigma^2, nrow = n), log=TRUE)
  ret <- sum(res)
  attr(ret,"components") <- res
  return(ret)
}

xthetalliknoODE <- function(x, CovV, CovR, sigma, y, grad = F)  {
  #a <- theta[1]
  #b <- theta[2]
  #c <- theta[3]
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  #f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  #frV <- (f[,1] - CovV$mphi %*% Vsm)
  
  # R
  #CovR <- calCov(phi[3:4])
  #frR <- (f[,2] - CovR$mphi %*% Rsm)
  
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2, 0,   - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2, 0,   - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  
  if(grad) {
    # V contrib
    # Vtemp <- diag( c*(1 - x[,1]^2)) - CovV$mphi
    # Rtemp <- diag( rep(c,n))
    # aTemp <- rep(0,n)
    # bTemp <- rep(0,n)
    # cTemp <- f[,1] / c
    # VC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovV$Kinv %*% frV
    
    # R contrib
    # Vtemp <- diag( rep( -1/c, n) )
    # Rtemp <- diag( rep( -b/c, n) ) - CovR$mphi
    # aTemp <- rep(1/c,n)
    # bTemp <- -Rsm/c
    # cTemp <- f[,2] * (-1/c)
    # RC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovR$Kinv %*% frR
    
    C3 <- c(2 * CovV$Cinv %*% Vsm,  2 * CovR$Cinv %*% Rsm)
    C1 <- c( 2 * (Vsm - y[,1]) / sigma^2 ,  2 * (Rsm - y[,2]) / sigma^2)
    
    #attr(ret,"grad") <- c((VC2 + RC2 + C3 + C1) * (-0.5))
    attr(ret,"grad") <- c((C3 + C1) * (-0.5))
  }
  
  return(ret)
  
  
}