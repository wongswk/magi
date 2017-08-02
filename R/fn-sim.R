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

getMeanDerivCurve <- function(x, y, dy, x.new, phi.mat, delta = 1e-9, sigma.mat, gamma.mat){
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
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    if(is.null(gamma.mat)){
      gamma <- 0
    }else{
      gamma <- gamma.mat[it]
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
    chol(C) %*% eta.mat[it,]
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
#### start of code ####

library(rstan)
fn.true <- read.csv("data/FN.csv")
fn.true$time <- seq(0,20,0.05)
matplot(fn.true$time, data.matrix(fn.true[,-3]), type="l", lty=1)

abc = c(0.2, 0.2, 3)

fn.true$dVtrue = with(fn.true, abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(fn.true, -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=0.1)


# fit_fit <- stan(file="stan/ode-fit.stan",
#                 data=list(N=nrow(fn.sim),
#                           robs=fn.sim$Rtrue,
#                           vobs=fn.sim$Vtrue),
#                 iter=200, chains=5)

fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=41),]
matplot(fn.sim$time, data.matrix(fn.sim[,-3]), type="l", lty=1)

gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           # drobs=fn.sim$dRtrue,
                           # dvobs=fn.sim$dVtrue,
                           time=fn.sim$time,
                           lambda=0.0083),
                 iter=100, chains=1)

traceplot(gpsmooth)
gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

plot(gpsmooth_ss$sigma, type="l",main="sigma")
abline(h=0.1, col=2)
hist(gpsmooth_ss$sigma, breaks = 50,main="sigma")
abline(v=0.1, col=2)

plot(gpsmooth_ss$gamma, type="l",main="gamma")
hist(gpsmooth_ss$gamma, breaks = 50,main="gamma")

vdRmcurve <- getMeanDerivCurve(x=fn.sim$time, y=fn.sim$Rtrue, dy=fn.sim$dRtrue, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

vdVmcurve <- getMeanDerivCurve(x=fn.sim$time, y=fn.sim$Vtrue, dy=fn.sim$dVtrue, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

Rpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$rphi, eta.mat = gpsmooth_ss$reta)
Vpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$vphi, eta.mat = gpsmooth_ss$veta)

dVdRpostsample <- getdVdR(abc.mat = gpsmooth_ss$abc, rtrue.mat = Rpostsample, vtrue.mat = Vpostsample)

matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.true$time, head(t(vdRmcurve),nrow(fn.true)), col="pink",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(Rpostsample), col="skyblue",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.true$time, head(t(vdRmcurve),nrow(fn.true)), col="pink",add=TRUE, type="l",lty=1)
matplot(fn.true$time, tail(t(vdRmcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(dVdRpostsample[,,"drobs"]), col="skyblue",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.true$time, head(t(vdVmcurve),nrow(fn.true)), col="pink",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(Vpostsample), col="skyblue",add=TRUE, type="l",lty=1)

matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.true$time, head(t(vdVmcurve),nrow(fn.true)), col="pink",add=TRUE, type="l",lty=1)
matplot(fn.true$time, tail(t(vdVmcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(dVdRpostsample[,,"dvobs"]), col="skyblue",add=TRUE, type="l",lty=1)

save(gpsmooth_ss, vdRmcurve, vdVmcurve, file="dump.RData")

gpsmooth_ss$abc
gpsmooth_ss$sigma
gpsmooth_ss$rphi
gpsmooth_ss$vphi

vdRmcurve[,seq(1,401,length=41)]
