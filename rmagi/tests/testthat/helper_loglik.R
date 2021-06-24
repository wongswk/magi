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
  n <- nrow(y)
  sigma <- tail(phisig,1)
  phiVR <- head(phisig, length(phisig)-1)
  res <- c(0,0)

  # V
  CovV <- calCov(head(phiVR, length(phiVR)/2), rInput, signrInput, kerneltype=kerneltype)
  Kv <- CovV$C+diag(sigma^2, nrow = n)
  Kv.l <- t(chol(Kv))
  Kv.l.inv <- solve(Kv.l)
  veta <- Kv.l.inv %*% y[,1]
  res[1] <- -n/2*log(2*pi) - sum(log(diag(Kv.l))) - 0.5*sum(veta^2)
  # R
  CovR <- calCov(tail(phiVR, length(phiVR)/2), rInput, signrInput, kerneltype=kerneltype)
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
    dVdphiAll <- apply(CovV$dCdphiCube, 3, function(dCdphiEach)
      sum(facVtemp*dCdphiEach)/2)

    # R contrib
    Kr.inv <- t(Kr.l.inv)%*%Kr.l.inv
    alphaR <- t(Kr.l.inv)%*%reta
    facRtemp <- alphaR%*%t(alphaR) - Kr.inv
    dRdsig <- sigma*sum(diag(facRtemp))
    dRdphiAll <- apply(CovR$dCdphiCube, 3, function(dCdphiEach)
      sum(facRtemp*dCdphiEach)/2)

    attr(ret,"grad") <- c(dVdphiAll, dRdphiAll, dVdsig+dRdsig)
  }
  return(ret)
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

#' marginal log likelihood for phi conditional on x and xdot
#'
#' used for fitting phi with ODE and conditioning on latent variables.
#' used to draw phi from Gibbs Sampler (possibly with HMC).
#'
#' this idea is not used in the end
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
