#### utility functions ####

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


