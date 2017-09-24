calCov <- function(phi, kerneltype="matern", rInput=NULL, signrInput=NULL) {
  if(is.null(rInput)){
    rInput <- r
  }
  if(is.null(signrInput)){
    signrInput <- signr
  }
  if(kerneltype=="matern"){
    ret <- calCov2(phi, rInput, signrInput)
  }else if(kerneltype=="rbf"){
    ret <- calCovRBF(phi, rInput, signrInput)
  }else if(kerneltype=="compact1"){
    ret <- calCovCompact1(phi, rInput, signrInput)
  }else{
    stop("kerneltype not specified correctly")
  }
  
  retmore <- with(ret, {
    C <- C + 1e-7 * diag( nrow(rInput))
    Cdecomp <- eigen(C)
    Ceigen1over <- 1/Cdecomp$value
    CeigenVec <- Cdecomp$vectors
    Cinv <- CeigenVec%*%(Ceigen1over*t(CeigenVec))
    mphi <-  Cprime %*% Cinv
    Kphi <- Cdoubleprime - (Cprime %*% Cinv %*% t(Cprime))  + 1e-7 * diag( nrow(rInput))
    Kdecomp <- eigen(Kphi)
    Keigen1over <- 1/Kdecomp$values
    KeigenVec <- Kdecomp$vectors
    Kinv <- KeigenVec%*%(Keigen1over*t(KeigenVec))
    mphiLeftHalf <- Cprime %*% CeigenVec
    list(Ceigen1over = Ceigen1over,
         CeigenVec = CeigenVec,
         Cinv = Cinv, 
         mphi = mphi, 
         Kphi = Kphi, 
         Keigen1over = Keigen1over,
         KeigenVec = KeigenVec,
         Kinv = Kinv,
         mphiLeftHalf = mphiLeftHalf)
  })
  c(ret, retmore)
}

calCov2 <- function(phi, r, signr, complexity=3) {
  r2 <- r^2
  C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
  if(complexity==0){
    return(list(C = C))
  }
  Cprime  <- (signr)* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
  Cdoubleprime <- (-phi[1] * (sqrt(5)/phi[2]) * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3))) + (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((10*sqrt(5)*r)/(3*phi[2]^3)))
  
  dCdphi <- list(
    C/phi[1],
    phi[1] * ( - ((sqrt(5)*r)/phi[2]^2) - ((10*r2)/(3*phi[2]^3))) * exp((-sqrt(5)*r)/phi[2]) + C * (sqrt(5)*r)/phi[2]^2
  )
  return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
}

calCovRBF <- function(phi, r, signr, complexity=3) {
  r2 <- r^2
  
  C <- phi[1] * exp(-r2/(2*phi[2]^2))
  if(complexity==0){
    return(list(C = C))
  }
  Cprime  <- signr * C * r / (phi[2]^2)
  Cdoubleprime <- C * (1/phi[2]^2 - r2 / phi[2]^4)
  dCdphi <- list(
    C/phi[1],
    C*r2/phi[2]^3
  )
  return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
}

calCovCompact1 <- function(phi, r, signr, complexity=3, D=3) {
  r2 <- r^2
  jsmooth <- floor(D/2)+2
  C <- phi[1] * pmax(1-(r/phi[2]),0)^(jsmooth+1) * ((jsmooth+1)*(r/phi[2])+1)
  
  if(complexity==0){
    return(list(C = C))
  }
  Cprime  <- signr * phi[1] * (jsmooth+1)/phi[2] * pmax(1-r/phi[2],0)^jsmooth * (jsmooth+2) * r/phi[2]
  Cdoubleprime <- phi[1] * pmax(1-r/phi[2],0)^(jsmooth-1) * (jsmooth+1) * (jsmooth+2) / phi[2]^2 * (1-r/phi[2]-r*jsmooth/phi[2])
  
  if(complexity==1){
    return(list(C = C, Cprime=Cprime, Cdoubleprime=Cdoubleprime))
  }
  
  dCdphi <- list(
    C/phi[1],
    phi[1] * pmax(1-r/phi[2],0)^jsmooth * r^2/phi[2]^3 * (jsmooth+1)*(jsmooth+2)
  )
  
  return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
}


fODE <- function(theta, x) {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  V <- x[,1]
  R <- x[,2]
  
  Vdt <- c * (V - V^3 / 3 + R)
  Rdt <- -1/c * ( V - a + b * R)
  
  return(cbind(Vdt, Rdt))
  #return(list(fit = res, Vfit = V, Rfit=R, Vfitdt = Vdt, Rfitdt = Rdt))
  
}

loglik <- function(x, theta, CovV, CovR, sigma, y, lambda=1)  {
  if(length(lambda) < 3) lambda <- c(lambda, rep(1, 3-length(lambda)))
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
  res[1,] <- c( (-0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma)) * lambda[1], 
                (-0.5 * as.numeric(determinant(CovV$Kphi)$modulus) - 0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr) * lambda[2],  
                (-0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm) * lambda[3])
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( (-0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma)) * lambda[1], 
                (-0.5 * as.numeric(determinant(CovR$Kphi)$modulus) - 0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr) * lambda[2],  
                (-0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm) * lambda[3])
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}

logliknoODE <- function(x, CovV, CovR, sigma, y)  {
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  #f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  #fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma), 0,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  #CovR <- calCov(phi[3:4])
  #fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma), 0,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}

library(mvtnorm)
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

xthetallik <- function(x, theta, CovV, CovR, sigma, y, grad = F, lambda = 1)  {
  if(length(lambda) < 3) lambda <- c(lambda, rep(1, 3-length(lambda)))
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
  
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 * lambda[1], 
                -0.5 * t(frV) %*% CovV$Kinv %*% frV * lambda[2],
                - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm * lambda[3])
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 * lambda[1],
                -0.5 * t(frR) %*% CovR$Kinv %*% frR * lambda[2],
                - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm * lambda[3])
  
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
    C1 <- c( 2 * (Vsm - y[,1]) / sigma^2 ,  2 * (Rsm - y[,2]) / sigma^2,0,0,0)
    
    attr(ret,"grad") <- c(((VC2 + RC2) * lambda[2] + C3 * lambda[3]+ C1 * lambda[1]) * (-0.5))
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

phisigllik <- function(phisig, y, grad = F, kerneltype="matern"){
  n <- nrow(y)
  sigma <- phisig[5]
  res <- c(0,0)
  
  # V 
  CovV <- calCov(phisig[1:2], kerneltype)
  Kv <- CovV$C+diag(sigma^2, nrow = n)
  Kv.l <- t(chol(Kv))
  Kv.l.inv <- solve(Kv.l)
  veta <- Kv.l.inv %*% y[,1]
  res[1] <- -n/2*log(2*pi) - sum(log(diag(Kv.l))) - 0.5*sum(veta^2)
  # R
  CovR <- calCov(phisig[3:4], kerneltype)
  Kr <- CovR$C+diag(sigma^2, nrow = n)
  Kr.l <- t(chol(Kr))
  Kr.l.inv <- solve(Kr.l)
  reta <- Kr.l.inv %*% y[,2]
  res[2] <- -n/2*log(2*pi) - sum(log(diag(Kr.l))) - 0.5*sum(reta^2)
  ret <- sum(res)
  attr(ret,"components") <- res
  
  if(grad) {
    # V contrib
    Kv.inv <- t(Kv.l.inv)%*%Kv.l.inv
    alphaV <- t(Kv.l.inv)%*%veta
    facVtemp <- alphaV%*%t(alphaV) - Kv.inv
    dVdsig <- sigma*sum(diag(facVtemp))
    dVdphi1 <- sum(facVtemp*CovV$dCdphi[[1]])/2
    dVdphi2 <- sum(facVtemp*CovV$dCdphi[[2]])/2
    
    # R contrib
    Kr.inv <- t(Kr.l.inv)%*%Kr.l.inv
    alphaR <- t(Kr.l.inv)%*%reta
    facRtemp <- alphaR%*%t(alphaR) - Kr.inv
    dRdsig <- sigma*sum(diag(facRtemp))
    dRdphi1 <- sum(facRtemp*CovR$dCdphi[[1]])/2
    dRdphi2 <- sum(facRtemp*CovR$dCdphi[[2]])/2
    
    attr(ret,"grad") <- c(dVdphi1, dVdphi2, dRdphi1, dRdphi2, dVdsig+dRdsig)
  }
  return(ret)
}

truncEigen <- function(eigenValues, frobeniusNormApprox = 0.99){
  frobeniusNorm <- sum(eigenValues^2)
  frobeniusNormTrunc <- cumsum(eigenValues^2)
  min(which(frobeniusNormTrunc/frobeniusNorm > frobeniusNormApprox))
}

truncCovByEigen <- function(gpCov, cKeep, kKeep){
  cKeepId <- (ncol(gpCov$CeigenVec)-cKeep+1):ncol(gpCov$CeigenVec)
  gpCov$Ceigen1over <- gpCov$Ceigen1over[cKeepId]
  gpCov$CeigenVec <- gpCov$CeigenVec[,cKeepId]

  kKeepId <- (ncol(gpCov$KeigenVec)-kKeep+1):ncol(gpCov$KeigenVec)  
  gpCov$Keigen1over <- gpCov$Keigen1over[kKeepId]
  gpCov$KeigenVec <- gpCov$KeigenVec[,kKeepId]
  
  gpCov
}
