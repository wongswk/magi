#' calculate stationary Gaussian process kernel
#' 
#' currently supports matern, rbf, compact1, 
#' also returns m_phi and other matrix which will be needed for ODE inference
#' 
#' @export
calCov <- function(phi, rInput, signrInput, bandsize = NULL, complexity=3, kerneltype="matern",
                   noiseInjection = 1e-7) {
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