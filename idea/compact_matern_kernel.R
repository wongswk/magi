calCov3 <- function(phi, r) {
  r2 <- r^2
  phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
}

r <- seq(0,5,0.01)
plot(r, calCov3(c(1,1), r), type="l")

compact_matern <- function(r, bandrange, nu) {
  pmax(1-r/bandrange, 0)^nu
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

x <- seq(-1,1,0.001)
r <- as.matrix(dist(x))
cov1 <- calCov3(c(1,1/4), r)
cov2 <- compact_matern(r,1,5/2)
cov3 <- cov1*cov2

y1 <- MASS::mvrnorm(10, mu=rep(0,length(x)), cov1)
matplot(x, t(y1), type="l")


y2 <- MASS::mvrnorm(10, mu=rep(0,length(x)), cov2)
matplot(x, t(y2), type="l")

x <- seq(0,5,0.05)
distm <- as.matrix(dist(x))
kdistm <- compact_matern(distm,1,5/2)
invkdistm <- solve(kdistm)
image(invkdistm)
plot.ts(invkdistm[,length(x)/2])
plot.ts(invkdistm[(length(x)*0.45):(length(x)*0.55),length(x)/2])

