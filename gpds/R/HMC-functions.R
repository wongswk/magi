#' calculate stationary Gaussian process kernel
#' 
#' currently supports matern, rbf, compact1, 
#' also returns m_phi and other matrix which will be needed for ODE inference
#' 
#' @export
calCov <- function(phi, rInput, signrInput, bandsize = NULL, complexity=3, kerneltype="matern",
                   noiseInjection = 1e-7, forceDiagKphi = FALSE) {
  if(kerneltype=="matern"){
    ret <- calCovMatern(phi, rInput, signrInput, complexity)
  }else if(kerneltype=="rbf"){
    ret <- calCovRBF(phi, rInput, signrInput, complexity)
  }else if(kerneltype=="compact1"){
    ret <- calCovCompact1(phi, rInput, signrInput, complexity)
  }else if(kerneltype=="periodicMatern"){
    ret <- calCovPeriodicWarpMatern(phi, rInput, signrInput, complexity)
  }else if(kerneltype=="generalMatern"){
    ret <- calCovGeneralMatern(phi, rInput, signrInput, complexity)
  }else if(kerneltype=="rationalQuadratic"){
    ret <- calCovRationalQuadratic(phi, rInput, signrInput, complexity)
  }else if(grepl("^generalMatern-([0-9.]+)$", kerneltype)){
    df = gsub("^generalMatern-([0-9.]+)$", "\\1", kerneltype)
    ret <- calCovGeneralMatern(phi, rInput, signrInput, complexity, df=as.numeric(df))
  }else{
    stop("kerneltype not specified correctly")
  }
  
  ret$mu <- rep(0, nrow(ret$C))
  ret$dotmu <- rep(0, nrow(ret$C))
  ret$C <- ret$C + noiseInjection * diag( nrow(rInput))
  
  if(complexity %in% c(0,1)){
    return(ret)
  }
  
  if(is.null(bandsize)){
    bandsize <- nrow(rInput)
  }
  
  retmore <- with(ret, {
    dCdphiCube <- sapply(dCdphi, identity, simplify = "array")
    Cdecomp <- eigen(C)
    Ceigen1over <- 1/Cdecomp$value
    CeigenVec <- Cdecomp$vectors
    Cinv <- CeigenVec%*%(Ceigen1over*t(CeigenVec))
    mphi <-  Cprime %*% Cinv
    Kright <- sqrt(Ceigen1over) * t(CeigenVec) %*% t(Cprime)
    Kphi <- Cdoubleprime - t(Kright)%*%Kright  + noiseInjection * diag( nrow(rInput))
    if(forceDiagKphi){
      Kphi <- mean(diag(Kphi)) * diag(nrow(Kphi))  
    }
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
         mphiLeftHalf = mphiLeftHalf,
         dCdphiCube = dCdphiCube)
  })
  retmore <- bandCov(retmore, bandsize)
  out <- c(ret, retmore)
  out$dCdphi <- NULL
  out
}

#' calculate Matern Gaussian process kernel
#' 
#' only calculate core part of C, Cprime, Cprimeprime, dCdphi etc.
#' 
#' @export
calCovMatern <- function(phi, r, signr, complexity=3) {
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

#' calculate general Matern Gaussian process kernel
#' 
#' only calculate core part of C, Cprime, Cprimeprime, dCdphi etc.
#' 
#' @export
calCovGeneralMatern <- function(phi, r, signr, complexity=3, df = 2.01) {
  r2 <- r^2
  
  x4bessel <- sqrt(2.0 * df) * r / phi[2];
  
  C <- phi[1] * 2^(1-df) * exp(-lgamma(df)) * x4bessel^df * besselK(x4bessel, df)
  C[r==0] <- phi[1]
  if(complexity==0){
    return(list(C = C))
  }
  
  dCdphi2 <- C * (df / x4bessel - (besselK(x4bessel, df-1) + besselK(x4bessel, df+1))/(2*besselK(x4bessel, df)))
  dCdphi2 <- dCdphi2 * (-sqrt(2.0 * df) * r / phi[2]^2)
  dCdphi2[is.na(dCdphi2)] <- 0
  
  
  Cprime  <- C * (df / x4bessel - (besselK(x4bessel, df-1) + besselK(x4bessel, df+1))/(2*besselK(x4bessel, df)))
  Cprime <- Cprime * sqrt(2.0 * df) * signr / phi[2]
  Cprime[is.na(Cprime)] <- 0
  Cprime <- -Cprime
  
  Cdoubleprime <- -phi[1] * 2^(1-df) * exp(-lgamma(df)) * 2.0 * df / phi[2]^2 * (
    df*(df-1)*x4bessel^(df-2)*besselK(x4bessel, df) - df*x4bessel^(df-1)*(besselK(x4bessel, df-1)+besselK(x4bessel, df+1))
    + x4bessel^df*(besselK(x4bessel, df-2)+2*besselK(x4bessel, df)+besselK(x4bessel, df+2))/4)
  diag(Cdoubleprime) <- phi[1] * 2^(1-df) * exp(-lgamma(df)) * 2.0 * df / phi[2]^2 * gamma(df-1) * 2^(df-2)
  # Cdoubleprime <- -Cdoubleprime
  
  dCdphi <- list(
    C/phi[1],
    dCdphi2
  )
  return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
}

#' calculate RBF Gaussian process kernel
#' 
#' only calculate core part of C, Cprime, Cprimeprime, dCdphi etc.
#' 
#' @export
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

#' calculate Compact1 Gaussian process kernel
#' 
#' only calculate core part of C, Cprime, Cprimeprime, dCdphi etc. 
#' See overleaf writeup for details
#' 
#' @export
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

#' calculate rational quadratic Gaussian process kernel
#' 
#' only calculate core part of C, Cprime, Cprimeprime, dCdphi etc.
#' 
#' @export
calCovRationalQuadratic <- function(phi, r, signr, complexity=3) {
  df = 0.01
  r2 <- r^2
  
  C <- phi[1] * (1 + r2/(2*phi[2]^2) / df)^(-df)
  if(complexity==0){
    return(list(C = C))
  }
  Cprime  <- -phi[1] * (-df) * (1 + r2/(2*phi[2]^2) / df)^(-df-1) * (r/ phi[2]^2 / df) * signr
  Cdoubleprime <- phi[1] * (-df) * (-df-1) * (1 + r2/(2*phi[2]^2) / df)^(-df-2) * (r/ phi[2]^2 / df)^2 +
    phi[1] * (-df) * (1 + r2/(2*phi[2]^2) / df)^(-df-1) * (1/ phi[2]^2 / df)
  Cdoubleprime <- -Cdoubleprime
  dCdphi <- list(
    C/phi[1],
    phi[1] * (-df) * (1 + r2/(2*phi[2]^2) / df)^(-df-1) * (-2 * r2/(2*phi[2]^3) / df)
  )
  return(list(C = C, Cprime = Cprime, Cdoubleprime = Cdoubleprime, dCdphi = dCdphi))
}

#' calculate linear Gaussian process kernel
#' 
#' @export
calCovLinear <- function(phi, x, complexity=3) {
  C <- phi[1] * (1 + outer(x, x)/phi[2]^2)
  
  if(complexity==0){
    return(list(C = C))
  }
}

#' calculate linear Gaussian process kernel
#' 
#' @export
calCovLinear <- function(phi, x, complexity=3) {
  C <- phi[1] * (1 + outer(x, x)/phi[2]^2)
  
  if(complexity==0){
    return(list(C = C))
  }
}

#' calculate Neural Network Gaussian process kernel
#' 
#' @export
calCovNeuralNetwork <- function(phi, x, complexity=3) {
  xtilde <- cbind(1, x)
  sigmaMat <- diag(c(phi[1], phi[2]^2))
  innerproduct <- xtilde%*%sigmaMat%*%t(xtilde)
  innerproductnormed <- diag(1/sqrt(1+2*diag(innerproduct)))%*%(2*innerproduct)%*%diag(1/sqrt(1+2*diag(innerproduct)))
  C <- 2/pi*asin(innerproductnormed)
  if(complexity==0){
    return(list(C = C))
  }
}

#' calculate modulated squared exponential Gaussian process kernel
#' 
#' @export
calCovModulatedRBF <- function(phi, x, complexity=3) {
  sigmaG <- phi[2]
  sigmaU <- phi[3]
  
  sigmaSqS <- 2*sigmaG^2 + sigmaG^4/sigmaU^2
  sigmaSqE <- 1/(2/sigmaG^2 + 1/sigmaU^2)
  sigmaSqM <- 2*sigmaU^2 + sigmaG^2
  
  C <- exp(-as.matrix(dist(x))^2/(2*sigmaSqS))
  C <- t(C * exp(-x^2/(2*sigmaSqM))) * exp(-x^2/(2*sigmaSqM))
  C <- phi[1] * C * (sqrt(sigmaSqE)/sigmaU)
  
  if(complexity==0){
    return(list(C = C))
  }
}

#' calculate Periodic Warpped Matern Gaussian process kernel
#' 
#' @export
calCovPeriodicWarpMatern <- function(phi, r, signr, complexity=0) {
  newr <- abs(sin(r*pi/phi[3])*2) # equivalent to cbind(sin(x*2*pi/phi[3]), cos(x*2*pi/phi[3]))
  maternCov <- calCovMatern(phi, newr, signr, complexity)
  maternCov$Cdoubleprime <- maternCov$Cdoubleprime * (cos(r*pi/phi[3])*2 * pi/phi[3])^2 -
    maternCov$Cprime * abs(sin(r*pi/phi[3])*2) * (pi/phi[3])^2 * signr
  maternCov$Cprime <- maternCov$Cprime * sign(sin(r*pi/phi[3])*2) * cos(r*pi/phi[3])*2 * pi/phi[3]
  maternCov$dCdphi[[3]] <- maternCov$C * sign(sin(r*pi/phi[3])*2) * cos(r*pi/phi[3])*2 * r*pi * -1/phi[3]^2
  maternCov
}

#' calculate fn model ODE derivatives
#' 
#' @export
fODE <- function(theta, x) {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  V <- x[,1]
  R <- x[,2]
  
  Vdt <- c * (V - V^3 / 3 + R)
  Rdt <- -1/c * ( V - a + b * R)
  
  return(cbind(Vdt, Rdt))
}

#' full log likelihood in R with contribution from each component
#' 
#' value is exact for x, theta, phi, sigma. we simply need to sample from this
#' if computation is feasible. phi info is contained in CovV, CovR
#' 
#' @export
loglik <- function(x, theta, CovV, CovR, sigma, y, lambda=1)  {
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
  res[1,] <- c( (-0.5 * sum((Vsm - y[,1])^2, na.rm = TRUE) / sigma^2 - n * log(sigma)) * lambda[1], 
                (-0.5 * as.numeric(determinant(CovV$Kphi)$modulus) - 0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr) * lambda[2],  
                (-0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm) * lambda[3])
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( (-0.5 * sum((Rsm - y[,2])^2, na.rm = TRUE) / sigma^2 - n * log(sigma)) * lambda[1], 
                (-0.5 * as.numeric(determinant(CovR$Kphi)$modulus) - 0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr) * lambda[2],  
                (-0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm) * lambda[3])
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}

#' full log likelihood in R with contribution from each component
#' 
#' value is exact for x, theta, phi, sigma. we simply need to sample from this
#' if computation is feasible. phi info is supplied explicitly
#' 
#' @export
loglikOrig <- function(x, theta, phi, sigma, y, rInput, signrInput, kerneltype = "matern")  {
  CovV <- calCov(phi[1:2], rInput, signrInput, complexity = 3, kerneltype=kerneltype)
  CovR <- calCov(phi[3:4], rInput, signrInput, complexity = 3, kerneltype=kerneltype)
  loglik(x, theta, CovV, CovR, sigma, y)
}

#' full log likelihood without ODE
#' 
#' used for Gaussian process smoothing.
#' full likelihood on latent x and phi sigma.
#' If x is marginalized out, we get phisigllik
#' 
#' @seealso phisigllik
#' 
#' @export
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

#' full log likelihood for latent x and theta in R
#' 
#' used for Gaussian process ODE inference on x and theta only. 
#' likelihood value is proportional to phi and sigma, so this cannot be used to
#' draw phi and sigma
#' 
#' @export
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
  
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2, na.rm = TRUE) / sigma^2 * lambda[1], 
                -0.5 * t(frV) %*% CovV$Kinv %*% frV * lambda[2],
                - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm * lambda[3])
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2, na.rm = TRUE) / sigma^2 * lambda[1],
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
    C1[is.na(C1)] <- 0
    
    attr(ret,"grad") <- c(((VC2 + RC2) * lambda[2] + C3 * lambda[3]+ C1 * lambda[1]) * (-0.5))
  }
  
  return(ret)
  
  
}

#' marginal log likelihood for phi and sigma for Gaussian process smoothing
#' 
#' used for Gaussian process smoothing without ODE. 
#' used to draw phi and sigma directly from marginal likelihood.
#' mathematically, this is the integration of logliknoODE
#' 
#' @export
phisigllik <- function(phisig, y, rInput, signrInput, grad = F, kerneltype="matern"){
  n <- NROW(y)
  y <- as.matrix(y)  
  
  sigma <- tail(phisig,1)
  phiVR <- matrix(head(phisig, length(phisig)-1), ncol=ncol(y))
  res <- rep(0, ncol(y))
  
  CovV <- list()
  Kv <- list()
  Kv.l <- list()
  Kv.l.inv <- list()
  veta <- list()
  
  for(j in 1:ncol(y)){
    # V 
    CovV[[j]] <- calCov(phiVR[,j], rInput, signrInput, kerneltype=kerneltype)
    Kv[[j]] <- CovV[[j]]$C+diag(sigma^2, nrow = n)
    Kv.l[[j]] <- t(chol(Kv[[j]]))
    Kv.l.inv[[j]] <- solve(Kv.l[[j]])
    veta[[j]] <- Kv.l.inv[[j]] %*% y[,j]
    res[j] <- -n/2*log(2*pi) - sum(log(diag(Kv.l[[j]]))) - 0.5*sum(veta[[j]]^2)
  }
  ret <- sum(res)
  attr(ret,"components") <- res
  
  if(grad) {
    dVdphiAll <- list()
    dVdsig <- list()
    for(j in 1:ncol(y)){
      # V contrib
      Kv.inv <- t(Kv.l.inv[[j]])%*%Kv.l.inv[[j]]
      alphaV <- t(Kv.l.inv[[j]])%*%veta[[j]]
      facVtemp <- alphaV%*%t(alphaV) - Kv.inv
      dVdsig[[j]] <- sigma*sum(diag(facVtemp))
      dVdphiAll[[j]] <- apply(CovV[[j]]$dCdphiCube, 3, function(dCdphiEach)
        sum(facVtemp*dCdphiEach)/2)
    }
    attr(ret,"grad") <- c(unlist(dVdphiAll), Reduce('+', dVdsig))
  }
  return(ret)
}

#' marginal log likelihood for phi conditional on x and xdot
#' 
#' used for fitting phi with ODE and conditioning on latent variables.
#' used to draw phi from Gibbs Sampler (possibly with HMC).
#' 
#' @export
phillikwithxdotx <- function(phi, x, xdot, rInput, signrInput, kerneltype="matern"){
  CovV <- calCov(phi, rInput, signrInput, kerneltype=kerneltype)
  bigCov <- rbind(cbind(CovV$C, t(CovV$Cprime)),
                  cbind(CovV$Cprime, CovV$Cdoubleprime))
  diag(bigCov) <- diag(bigCov) + 1e-7
  
  bigCov.l <- t(chol(bigCov))
  eta <- solve(bigCov.l, c(x, xdot))
  
  ret <- -length(c(x, xdot))/2*log(2*pi) - sum(log(diag(bigCov.l))) - 0.5*sum(eta^2)
  
  return(ret)
}

#' calculate number of eigen values to preserve based on frobenius norm
#' @export
truncEigen <- function(eigenValues, frobeniusNormApprox = 0.99){
  frobeniusNorm <- sum(eigenValues^2)
  frobeniusNormTrunc <- cumsum(eigenValues^2)
  min(which(frobeniusNormTrunc/frobeniusNorm > frobeniusNormApprox))
}

#' truncate gpCov object for low rank approximation
#' 
#' the largest 1-over-eigenvalue will be preserved, and the rest will be deleted
#' for a low-rank representation from spectral decomposition
#' 
#' mphi SVD not helping because complexity is 2mn, comparing to original n^2
#' however we need m to be around n/2 to preserve accuracy due to low decrease d
#' 
#' @param cKeep number of eigen values to keep for C matrix
#' @param kKeep number of eigen values to keep for K matrix
#' @export
truncCovByEigen <- function(gpCov, cKeep, kKeep){
  cKeepId <- (ncol(gpCov$CeigenVec)-cKeep+1):ncol(gpCov$CeigenVec)
  gpCov$Ceigen1over <- gpCov$Ceigen1over[cKeepId]
  gpCov$CeigenVec <- gpCov$CeigenVec[,cKeepId]

  kKeepId <- (ncol(gpCov$KeigenVec)-kKeep+1):ncol(gpCov$KeigenVec)  
  gpCov$Keigen1over <- gpCov$Keigen1over[kKeepId]
  gpCov$KeigenVec <- gpCov$KeigenVec[,kKeepId]
  
  # mKeepId <- 1:mKeep
  # gpCov$mphiu <- gpCov$mphiu[mKeepId,]
  # gpCov$mphid <- gpCov$mphid[mKeepId]
  # gpCov$mphiv <- gpCov$mphiv[mKeepId,]
  
  gpCov
}

#' calculate number of eigen values to preserve based on frobenius norm
#' @export
bandCov <- function(gpCov, bandsize = 20){
  gpCov$CinvBand <- mat2band(gpCov$Cinv, bandsize)
  gpCov$mphiBand <- mat2band(gpCov$mphi, bandsize)
  gpCov$KinvBand <- mat2band(gpCov$Kinv, bandsize)
  gpCov$bandsize <- bandsize
  gpCov
}

#' Converts a matrix to banded form
#' 
#' @export
mat2band <- function(a, bandsize) {
  N <- nrow(a)
  A <- matrix(0,nrow = 2*bandsize+1, ncol = N)
  
  for (j in 1:N) {
    k <- bandsize + 1 - j
    for (i in max(1,j-bandsize):min(N,j+bandsize)) {
      A[k+i,j] <- a[i,j] 
    }
  }
  
  A
}


#' get posterior mean curve for value conditioning on observed y, phi, sigma
#' 
#' use to visulize Gaussian process smoothing without ODE, and for initialize 
#' the x theta sampler
#' 
#' @export
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

#' get mphi in the same format of calCov return values
#' 
#' @export
getCovMphi <- function(kernel, xsim, xsim.obs, config=list()){
  bandsize <- config$bandsize
  forceDiagKphi <- config$forceDiagKphi
  if(is.null(bandsize)){
    bandsize <- nrow(xsim)
  }
  if(is.null(forceDiagKphi)){
    forceDiagKphi <- FALSE
  }
  if(kernel=="finiteDifference2h"){
    diffMat <- matrix(0, nrow=nrow(xsim), ncol=nrow(xsim))
    for(i in 1:nrow(xsim)){
      if(i==1){
        diffMat[i,1] <- -1
        diffMat[i,2] <- 1
        diffMat[i,] <- diffMat[i,] / (xsim$time[2] - xsim$time[1])
      }else if(i==nrow(xsim)){
        diffMat[i,nrow(xsim)] <- 1
        diffMat[i,nrow(xsim)-1] <- -1
        diffMat[i,] <- diffMat[i,] / (xsim$time[nrow(xsim)] - xsim$time[nrow(xsim)-1])
      }else{
        h1 <- xsim$time[i+1] - xsim$time[i]
        h2 <- xsim$time[i] - xsim$time[i-1]
        diffMat[i,i-1] <- -h1^2
        diffMat[i,i] <- h1^2 - h2^2
        diffMat[i,i+1] <- h2^2
        diffMat[i,] <- diffMat[i,] / (h1*h2^2 + h1^2*h2)
      }
    }
    curCovNew <- rep(list(list()), ncol(xsim)-1)
    for(j in 1:length(curCovNew)){
      curCovNew[[j]]$mphi <- diffMat
      curCovNew[[j]]$mphiBand <- mat2band(diffMat, bandsize=bandsize)
      curCovNew[[j]]$Kinv <- diag(nrow(xsim))
      curCovNew[[j]]$KinvBand <- mat2band(curCovNew[[j]]$Kinv, bandsize=bandsize)
      curCovNew[[j]]$Cinv <- matrix(0, nrow(xsim), nrow(xsim))
      curCovNew[[j]]$CinvBand <- mat2band(curCovNew[[j]]$Cinv, bandsize=bandsize)
      curCovNew[[j]]$bandsize <- bandsize
    }
    return(curCovNew)
  }else if(kernel=="finiteDifference1h"){
    diffMat <- matrix(0, nrow=nrow(xsim), ncol=nrow(xsim))
    for(i in 1:nrow(xsim)){
      if(i==1){
        diffMat[i,1] <- -1
        diffMat[i,2] <- 1
        diffMat[i,] <- diffMat[i,] / (xsim$time[2] - xsim$time[1])
      }else if(i==nrow(xsim)){
        diffMat[i,nrow(xsim)] <- 1
        diffMat[i,nrow(xsim)-1] <- -1
        diffMat[i,] <- diffMat[i,] / (xsim$time[nrow(xsim)] - xsim$time[nrow(xsim)-1])
      }else{
        diffMat[i,i] <- -1
        diffMat[i,i+1] <- 1
        diffMat[i,] <- diffMat[i,] / (xsim$time[i+1] - xsim$time[i])
      }
    }
    curCovNew <- rep(list(list()), ncol(xsim)-1)
    for(j in 1:length(curCovNew)){
      curCovNew[[j]]$mphi <- diffMat
      curCovNew[[j]]$mphiBand <- mat2band(diffMat, bandsize=bandsize)
      curCovNew[[j]]$Kinv <- diag(nrow(xsim))
      curCovNew[[j]]$KinvBand <- mat2band(curCovNew[[j]]$Kinv, bandsize=bandsize)
      curCovNew[[j]]$Cinv <- matrix(0, nrow(xsim), nrow(xsim))
      curCovNew[[j]]$CinvBand <- mat2band(curCovNew[[j]]$Cinv, bandsize=bandsize)
      curCovNew[[j]]$bandsize <- bandsize
    }
    return(curCovNew)
  }else if(kernel=="zero"){
    curCovNew <- rep(list(list()), 2)
    for(j in 1:length(curCovNew)){
      curCovNew[[j]]$mphi <- matrix(0, nrow(xsim), nrow(xsim))
      curCovNew[[j]]$mphiBand <- mat2band(matrix(0, nrow(xsim), nrow(xsim)), bandsize=bandsize)
      curCovNew[[j]]$mphiLeftHalf <- NULL
    }
    return(curCovNew)
  }
  
  foo.nobs <- outer(xsim.obs$time, t(xsim.obs$time),'-')[,1,]
  r.nobs <- abs(foo.nobs)
  signr.nobs <- -sign(foo.nobs)
  
  cursigma <- rep(NA, ncol(xsim)-1)
  curphi <- matrix(NA, 2, ncol(xsim)-1)
  
  for(j in 1:(ncol(xsim)-1)){
    priorFactor <- getFrequencyBasedPrior(xsim.obs[,1+j])
    
    desiredMode <- priorFactor["meanFactor"]
    betaRate <- uniroot(function(betaRate) pgamma(1, 1 + desiredMode*betaRate, betaRate)-0.95,
                        c(1e-3, 1e3))$root
    alphaRate <- 1 + desiredMode*betaRate
    
    fn <- function(par) {
      marlik <- phisigllik( par, data.matrix(xsim.obs[,1+j]), r.nobs, signr.nobs, grad=FALSE, kerneltype = kernel)
      penalty <- dnorm(par[2], max(xsim.obs$time)*priorFactor["meanFactor"], 
                       max(xsim.obs$time)*priorFactor["sdFactor"], log=TRUE)
      -(as.numeric(marlik) + penalty)
    }
    gr <- function(par) {
      marlik <- phisigllik( par, data.matrix(xsim.obs[,1+j]), r.nobs, signr.nobs, grad=TRUE, kerneltype = kernel)
      grad <- -as.vector(attr(marlik, "grad"))
      penalty <- (par[2] - max(xsim.obs$time)*priorFactor["meanFactor"]) / (max(xsim.obs$time)*priorFactor["sdFactor"])^2
      grad[2] <- grad[2] + penalty
      grad
    }
    testthat::expect_equal(gr(c(5,50,1))[2], (fn(c(5,50+1e-6,1)) - fn(c(5,50,1)))/1e-6, tolerance=1e-3)
    marlikmap <- optim(c(sd(xsim.obs[,1+j])/2, max(xsim.obs$time)/2, sd(xsim.obs[,1+j])/2), 
                       fn, gr, method="L-BFGS-B", lower = 0.0001,
                       upper = c(Inf, Inf, Inf))
    
    cursigma[j] <- marlikmap$par[3]
    curphi[,j] <- marlikmap$par[1:2]
  }
  
  foo <- outer(xsim$time, t(xsim$time),'-')[,1,]
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], abs(foo), -sign(foo), kerneltype=kernel, 
                      bandsize = bandsize, forceDiagKphi = forceDiagKphi)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  curCov$cursigma <- cursigma
  curCov$curphi <- curphi
  curCov
}
