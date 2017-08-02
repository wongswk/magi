
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

##   phiV is phi[1:2], phiR is phi[3:4]
calCov <- function(phi) {
  C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
  Cprime  <- (signr)* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
  Cdoubleprime <- (-phi[1] * (sqrt(5)/phi[2]) * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3))) + (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((10*sqrt(5)*r)/(3*phi[2]^3)))
  
  C <- C + 1e-9 * diag( nrow(r))
  
  Cinv <- solve(C)
  mphi <-  Cprime %*% Cinv
  Kphi <- Cdoubleprime - (Cprime %*% Cinv %*% t(Cprime))  + 1e-9 * diag( nrow(r))
  
  return(list(C = C, Cinv = Cinv, mphi = mphi, Kphi = Kphi))
}

loglik <- function(x, theta, phi, sigma, y)  {
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
  CovV <- calCov(phi[1:2])
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovV$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovR$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}


# Set up for 41 points
tvec41 <- seq(0,20, by = 0.05*10)  
foo <- outer(tvec41, t(tvec41),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)


# Calc log-lik of STAN output posterior means
loglik(cbind(colMeans(vdVmcurve[,seq(1,401,length=41)]),colMeans(vdRmcurve[,seq(1,401,length=41)])), colMeans(gpsmooth_ss$abc), c(colMeans(gpsmooth_ss$vphi), colMeans(gpsmooth_ss$rphi)),mean(gpsmooth_ss$sigma), fn.sim[,1:2])

# "best" log-likelihood based on truth. phi vector found by optim with other inputs set at truth
loglik( VRtrue[seq(1,401,length=41),], c(0.2,0.2,3), c(1.9840824, 1.1185157, 0.9486433, 3.2682434), 0.1,  fn.sim[,1:2])
