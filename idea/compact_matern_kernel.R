calCov3 <- function(phi, r) {
  r2 <- r^2
  phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
}

r <- seq(0,5,0.01)
plot(r, calCov3(c(1,1), r), type="l")

compact_matern <- function(r, bandrange, nu) {
  pmax(1-r/bandrange, 0)^nu
}

compact_smooth <- function(r, D, q) {
  j <- floor(D/2)+q+1
  if(q==0){
    ret <- pmax(1-r,0)^j
  }else if(q==1){
    ret <- pmax(1-r,0)^(j+1) * ((j+1)*r+1)
  }else if(q==2){
    ret <- pmax(1-r,0)^(j+2) * ((j^2+4*j+3)*r^2+(3*j+6)*r+3)/3
  }else if(q==3){
    ret <- pmax(1-r,0)^(j+3) * ((j^3+9*j^2+23*j+15)*r^3+(6*j^2+36*j+45)*r^2+(15*j+45)*r+15)/15
  }else{
    stop("invalid q")
  }
  return(ret)
}

r <- seq(0,5,0.01)
plot(r, compact_matern(r,1,5/2), type="l")

plot(r, compact_matern(r,5,5/2)*calCov3(c(1,5), r), type="l")

r <- seq(0,1,0.001)
plot(r, compact_matern(r,1,5/2)*calCov3(c(1,1), r), type="l")
lines(r, compact_matern(r,1,5/2), col=2)
lines(r, calCov3(c(1,1), r), col=3)
lines(r, calCov3(c(1,1/3), r), col=4)
lines(r, calCov3(c(1,1/4), r), col=5)
lines(r, compact_smooth(r,3,1), col=6)

x <- seq(-1,1,0.01)
r <- as.matrix(dist(x))
cov1 <- calCov3(c(1,1/4), r)
cov2 <- compact_matern(r,1,5/2)

y1 <- MASS::mvrnorm(10, mu=rep(0,length(x)), cov1)
matplot(x, t(y1), type="l")


y2 <- MASS::mvrnorm(10, mu=rep(0,length(x)), cov2)
matplot(x, t(y2), type="l")

cov3 <- compact_smooth(r,3,1)
y3 <- MASS::mvrnorm(10, mu=rep(0,length(x)), cov3)
matplot(x, t(y3), type="l")

x <- seq(0,5,0.05)
distm <- as.matrix(dist(x))
kdistm <- compact_smooth(distm,3,1)
min(eigen(kdistm)$value)
invkdistm <- solve(kdistm)
plot(x, kdistm[,length(x)/2], type="l")
plot(x, invkdistm[,length(x)/2], type="l")
plot.ts(invkdistm[(length(x)*0.45):(length(x)*0.55),length(x)/2])

