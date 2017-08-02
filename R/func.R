#### utility functions ####
getMeanCurve <- function(x, y, x.new, phi.mat, delta = 1e-9, sigma.mat){
  tvec <- c(x.new,x)
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  
  signr <- -sign(foo)
  
  t(sapply(1:nrow(phi.mat), function(it){
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    diag(C) <- diag(C)+delta
    diag(C)[-(1:length(x.new))] <- diag(C)[-(1:length(x.new))]+sigma^2
    C[1:length(x.new),-(1:length(x.new))]%*%solve(C[-(1:length(x.new)),-(1:length(x.new))], y)
  }))
}

getDerivCurve <- function(x, y, dy, x.new, phi.mat, delta = 1e-9, sigma.mat, gamma.mat){
  tvec.d <- c(x.new,x)
  tvec.v <- x
  tvec <- c(tvec.d, tvec.v)
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  
  signr <- -sign(foo)
  
  t(sapply(1:nrow(phi.mat), function(it){
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    gamma <- gamma.mat[it]
    
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    diag(C) <- diag(C)
    C <- C[-(1:length(tvec.d)),-(1:length(tvec.d))]
    
    Cprime  <- signr* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
    Cprime <- Cprime[(1:length(tvec.d)),-(1:length(tvec.d))]
    Cdoubleprime <- (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((5*sqrt(5)*r)/(3*phi[2]^3)) - ((25*r2)/(3*phi[2]^4)))
    Cdoubleprime <- Cdoubleprime[(1:length(tvec.d)),(1:length(tvec.d))]
    Cdoubleprime[-(1:length(x.new))] <- Cdoubleprime[-(1:length(x.new))] + gamma
    
    M <- rbind(cbind(Cdoubleprime, Cprime),
               cbind(t(Cprime), C+sigma))
    
    M[1:length(x.new),-(1:length(x.new))]%*%solve(M[-(1:length(x.new)),-(1:length(x.new))], c(dy, y))
  }))
}

getDerivCurve2 <- function(x, y, x.new, phi.mat, delta = 1e-9, sigma.mat, gamma.mat){
  tvec.d <- c(x.new)
  tvec.v <- x
  tvec <- c(tvec.d, tvec.v)
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  
  signr <- -sign(foo)
  
  t(sapply(1:nrow(phi.mat), function(it){
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    gamma <- gamma.mat[it]
    
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    C <- C[-(1:length(tvec.d)),-(1:length(tvec.d))]
    
    Cprime  <- signr* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
    Cprime <- Cprime[(1:length(tvec.d)),-(1:length(tvec.d))]
    
    Cdoubleprime <- (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((5*sqrt(5)*r)/(3*phi[2]^3)) - ((25*r2)/(3*phi[2]^4)))
    Cdoubleprime <- Cdoubleprime[(1:length(tvec.d)),(1:length(tvec.d))]
    Cdoubleprime[-(1:length(x.new))] <- Cdoubleprime[-(1:length(x.new))] + gamma^2
    
    M <- rbind(cbind(Cdoubleprime, Cprime),
               cbind(t(Cprime), C+sigma^2))
    
    M[1:length(x.new),-(1:length(x.new))]%*%solve(M[-(1:length(x.new)),-(1:length(x.new))], y)
  }))
}

getMeanDerivCurve <- function(x, y.mat, dy.mat, x.new, phi.mat, delta = 1e-9, sigma.mat, gamma.mat){
  tvec <- c(x.new,x,x.new,x)
  id.dnew <- 1:length(x.new)
  id.dobs <- 1:length(x) + length(x.new)
  id.vnew <- length(tvec)/2 + 1:length(x.new)
  id.vobs <- tail(1:length(tvec), length(x))
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  
  signr <- -sign(foo)
  
  t(sapply(1:nrow(phi.mat), function(it){
    y <- y.mat[it,]
    dy <- dy.mat[it,]
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    if(is.null(gamma.mat)){
      gamma <- 0
    }else{
      gamma <- gamma.mat[it]
    }
    
    if(is.null(sigma.mat)){
      sigma <- 0.1
    }else{
      sigma <- sigma.mat[it]
    }
    
    
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    Cprime  <- signr* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
    Cdoubleprime <- (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((5*sqrt(5)*r)/(3*phi[2]^3)) - ((25*r2)/(3*phi[2]^4)))
    
    M <- matrix(NA, ncol=length(tvec), nrow=length(tvec))
    
    M[c(id.dnew,id.dobs),c(id.dnew,id.dobs)] <- Cdoubleprime[c(id.dnew,id.dobs),c(id.dnew,id.dobs)]
    M[c(id.vnew,id.vobs),c(id.vnew,id.vobs)] <- C[c(id.vnew,id.vobs),c(id.vnew,id.vobs)]
    
    M[c(id.dnew,id.dobs),c(id.vnew,id.vobs)] <- Cprime[c(id.dnew,id.dobs),c(id.vnew,id.vobs)]
    M[c(id.vnew,id.vobs),c(id.dnew,id.dobs)] <- t(Cprime[c(id.dnew,id.dobs),c(id.vnew,id.vobs)])
    
    diag(M)[id.vobs] <- diag(M)[id.vobs]+sigma^2
    diag(M)[id.dobs] <- diag(M)[id.dobs]+gamma^2
    
    stopifnot(is.finite(M))
    diag(M) <- diag(M)+delta
    
    M[c(id.vnew,id.dnew),c(id.vobs,id.dobs)]%*%solve(M[c(id.vobs,id.dobs),c(id.vobs,id.dobs)], c(y,dy))
  }))
}

getX <- function(r, phi.mat, eta.mat, delta = 1e-9){
  r2 <- r^2
  t(sapply(1:nrow(phi.mat), function(it){
    phi <- phi.mat[it,]
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    diag(C) <- diag(C) + delta
    t(chol(C)) %*% eta.mat[it,]
  }))
}

getdVdR <- function(abc.mat, rtrue.mat, vtrue.mat){
  result <- sapply(1:nrow(abc.mat), function(it){
    abc <- abc.mat[it,]
    vtrue <- vtrue.mat[it,]
    rtrue <- rtrue.mat[it,]
    
    dvobs = abc[3] * (vtrue - (vtrue^3)/3.0 + rtrue)
    drobs = -1.0/abc[3] * (vtrue - abc[1] + abc[2]*rtrue)
    cbind(drobs, dvobs)
  }, simplify = "array")
  aperm(result, c(3,1,2))
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

##   phiV is phi[1:2], phiR is phi[3:4]
calCov <- function(phi, r) {
  r2 <- r^2
  C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
  Cprime  <- (signr)* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
  Cdoubleprime <- (-phi[1] * (sqrt(5)/phi[2]) * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3))) + (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((10*sqrt(5)*r)/(3*phi[2]^3)))
  
  C <- C + 1e-9 * diag( nrow(r))
  
  Cinv <- solve(C)
  mphi <- Cprime %*% Cinv
  Kphi <- Cdoubleprime - (Cprime %*% Cinv %*% t(Cprime))  + 1e-9 * diag( nrow(r))
  
  return(list(C = C, Cinv = Cinv, mphi = mphi, Kphi = Kphi, Cprime=Cprime, Cdoubleprime=Cdoubleprime))
}

loglik <- function(x, theta, phi, sigma, y, r)  {
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
  CovV <- calCov(phi[1:2], r)
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovV$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  CovR <- calCov(phi[3:4], r)
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovR$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}
