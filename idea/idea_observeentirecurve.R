sourceCpp("../src/wrapper.cpp")
byl = 0.01
x <- seq(0,1,byl)
xdist <- as.matrix(dist(x))
yobs <- cbind(sin(x), cos(x))

#' without error term, the likelihood is numerically every unstable,
#' use error 0.01 for this example
fn <- function(par) -phisigllikTest( c(par,1e-2), yobs, xdist)$value
gr <- function(par) -as.vector(phisigllikTest( c(par,1e-2), yobs, xdist)$grad)[1:4]
x0 <- rep(1,4)
marlikmap <- optim(rep(1,4), fn, gr, method="L-BFGS-B", lower = 1e-10)

marlikmap$par
marlikmap$value

for(i in 1:4){
  x1 = x0
  x1[i] = x1[i] + 1e-5
  print((fn(x1) - fn(x0))/1e-5)
}
gr(x0)



thetasigma <- marlikmap$par
gr(thetasigma)
thetasigma[5] <- 0.01

thetasigma1 <- thetasigma
thetasigma2 <- thetasigma
thetasigma2 <- thetasigma1*1.1

byl = 0.01/16
x <- seq(0,1,byl)
xdist <- as.matrix(dist(x))
yobs <- cbind(sin(x), cos(x))


phisigllikTest( thetasigma1, yobs, xdist)$value - phisigllikTest( thetasigma2, yobs, xdist)$value
# log likelihood difference explode


#### brownian motion mean parameter ####
#' with inf observations the parameter cannot be pin down exactly
rm(list=ls())
sigmasq <- 1
time <- seq(0,1,length=2^20+1)
mu <- 1
ytrue <- rnorm(length(time)-1, diff(time)*mu, sqrt(diff(time)*sigmasq))



sum(dnorm(ytrue, diff(time)*1.0, sqrt(diff(time)*sigmasq), log=TRUE)) - 
  sum(dnorm(ytrue, diff(time)*0.9, sqrt(diff(time)*sigmasq), log=TRUE))

time.obs <- time
y.obs <- ytrue
for(coarse in seq(20,1)){
  time.obs <- time.obs[seq(1,length(time.obs),2)]
  y.obs <- y.obs[seq(1,length(y.obs),2)] + y.obs[seq(1,length(y.obs),2)+1]
  print(sum(dnorm(y.obs, diff(time.obs)*1.0, sqrt(diff(time.obs)*sigmasq), log=TRUE)) - 
          sum(dnorm(y.obs, diff(time.obs)*0.9, sqrt(diff(time.obs)*sigmasq), log=TRUE)))
}
