#' Calculate stationary Gaussian process kernel
#' 
#' Covariance calculations for Gaussian process kernels.
#' Currently supports matern, rbf, compact1, periodicMatern, generalMatern, and rationalQuadratic kernels.
#' Can also return m_phi and other additional quantities useful for ODE inference.
#'
#' @param phi the kernel hyper-parameters. See details for hyper-parameter specification for each \code{kerneltype}.
#' @param rInput the distance matrix between all time points s and t, i.e., |s - t|
#' @param signrInput the sign matrix of the time differences, i.e., sign(s - t)
#' @param bandsize size for band matrix approximation. See details.
#' @param complexity integer value for the complexity of the kernel calculations desired:
#' \itemize{
#'   \item 0 includes C only
#'   \item 1 additionally includes Cprime, Cdoubleprime, dCdphi
#'   \item 2 or above additionally includes Ceigen1over, CeigenVec, Cinv, mphi, Kphi, Keigen1over, KeigenVec, Kinv, mphiLeftHalf, dCdphiCube
#' }
#' See details for their definitions.
#' @param kerneltype must be one of \code{matern}, \code{rbf}, \code{compact1}, \code{periodicMatern}, \code{generalMatern}, \code{rationalQuadratic}. See details for the kernel formulae.
#' @param df degrees of freedom, for \code{generalMatern} and \code{rationalQuadratic} kernels only.  Default is \code{df=2.01} for \code{generalMatern} and \code{df=0.01} for \code{rationalQuadratic}.
#' @param noiseInjection a small value added to the diagonal elements of C and Kphi for numerical stability
#' 
#' @return A list containing the kernel calculations included by the value of \code{complexity}.
#' 
#' @details 
#' The covariance formulae and the hyper-parameters \code{phi} for the supported kernels are as follows.  Stationary kernels have \eqn{C(s,t) = C(r)} where \eqn{r = |s-t|} is the distance between the two time points.   Generally, the hyper-parameter \code{phi[1]} controls the overall variance level while \code{phi[2]} controls the bandwidth.
#' \describe{
#'   \item{\code{matern}}{ This is the simplified Matern covariance with \code{df = 5/2}:
#'   \deqn{C(r) = phi[1] * (1 + \sqrt 5 r/phi[2] + 5r^2/(3 phi[2]^2)) * \exp(-\sqrt 5 r/phi[2])}
#'   }
#'   \item{\code{rbf}}{
#'   \deqn{C(r) = phi[1] * \exp(-r^2/(2 phi[2]^2))}
#'   }
#'   \item{\code{compact1}}{
#'   \deqn{C(r) = phi[1] * \max(1-r/phi[2],0)^4 * (4r/phi[2]+1) }
#'   }
#'   \item{\code{periodicMatern}}{
#'   Define \eqn{r' =  | \sin(r \pi/phi[3])*2 |}.  Then the covariance is given by \eqn{C(r')} using the Matern formula.
#'   }
#'   \item{\code{generalMatern}}{
#'   \deqn{C(r) = phi[1] * 2^(1-df) / \Gamma(df) * ( \sqrt(2.0 * df) * r / phi[2] )^df * besselK( \sqrt(2.0 * df) * r / phi[2] , df)}
#'   where \code{besselK} is the modified Bessel function of the second kind.
#'   }
#'   \item{\code{rationalQuadratic}}{
#'   \deqn{C(r) = phi[1] * (1 + r^2/(2 df phi[2]^2))^(-df)}
#'   }
#' }
#' 
#' The kernel calculations available and their definitions are as follows: 
#' \describe{
#'   \item{C}{The covariance matrix corresponding to the distance matrix \code{rInput}.}
#'   \item{Cprime}{The cross-covariance matrix  \eqn{d C(s,t) / ds}.}
#'   \item{Cdoubleprime}{The cross-covariance matrix  \eqn{d^2 C(s,t) / ds dt}.}
#'   \item{dCdphi}{A list with the matrices \eqn{dC / dphi} for each element of phi.}
#'   \item{Ceigen1over}{The reciprocals of the eigenvalues of C.}
#'   \item{CeigenVec}{Matrix of eigenvectors of C.}
#'   \item{Cinv}{The inverse of C.}
#'   \item{mphi}{The matrix \code{Cprime * Cinv}.}
#'   \item{Kphi}{The matrix \code{Cdoubleprime - Cprime * Kinv * t(Cprime)}.}
#'   \item{Keigen1over}{The reciprocals of the eigenvalues of Kphi.}
#'   \item{Kinv}{The inverse of Kphi.}
#'   \item{mphiLeftHalf}{The matrix \code{Cprime * CeigenVec}.}
#'   \item{dCdphiCube}{\eqn{dC / dphi} as a 3-D array, with the third dimension corresponding to the elements of phi.}
#' }
#' 
#' If \code{bandsize} is a positive integer, additionally CinvBand, mphiBand, and KinvBand are provided in the return list, which are
#' band matrix approximations to Cinv, mphi, and Kinv with the specified \code{bandsize}.
#' 
#' 
#' @examples 
#' foo  <- outer(0:40, t(0:40),'-')[,1,]
#' r <- abs(foo)
#' signr <- -sign(foo)
#' calCov(c(0.2, 2), r, signr, bandsize=20, kerneltype="generalMatern", df=2.01)
#' 
#' @export
calCov <- function(phi, rInput, signrInput, bandsize = NULL, complexity=3, kerneltype="matern", df,
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
    if (missing(df))
      ret <- calCovGeneralMatern(phi, rInput, signrInput, complexity)
    else
      ret <- calCovGeneralMatern(phi, rInput, signrInput, complexity, df)
  }else if(kerneltype=="rationalQuadratic"){
    if (missing(df))
      ret <- calCovRationalQuadratic(phi, rInput, signrInput, complexity)
    else
      ret <- calCovRationalQuadratic(phi, rInput, signrInput, complexity, df)
  }else{
    stop("kerneltype not specified correctly")
  }
  
  ret$mu <- rep(0, nrow(ret$C))
  ret$dotmu <- rep(0, nrow(ret$C))
  ret$C <- ret$C + noiseInjection * diag( nrow(rInput))
  ret$tvecCovInput <- rInput[1,]

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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
calCovRationalQuadratic <- function(phi, r, signr, complexity=3, df=0.01) {
  #df = 0.01
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
#' @noRd
calCovLinear <- function(phi, x, complexity=3) {
  C <- phi[1] * (1 + outer(x, x)/phi[2]^2)
  
  if(complexity==0){
    return(list(C = C))
  }
}

#' calculate Neural Network Gaussian process kernel
#'
#' @noRd
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
#' @noRd
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
#' @noRd
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
#'
#' @noRd
bandCov <- function(gpCov, bandsize = 20){
  gpCov$CinvBand <- mat2band(gpCov$Cinv, bandsize)
  gpCov$mphiBand <- mat2band(gpCov$mphi, bandsize)
  gpCov$KinvBand <- mat2band(gpCov$Kinv, bandsize)
  gpCov$bandsize <- bandsize
  gpCov
}

#' Converts a matrix to banded form
#'
#' @noRd
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
