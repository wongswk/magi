x <- diag(100)
diag(x[-1,-10]) <- 0.5
diag(x[-10,-1]) <- 0.5
xinv <- solve(x)
plot(abs(xinv[,1]))
plot(abs(xinv[1:20,1])) # inverse is compact as well
abs(xinv[-(1:10),1])

rm(list=ls())

source("../R/HMC-functions.R")
x <- seq(0,5,0.05)
distm <- as.matrix(dist(x))
kdistm <- calCovCompact1(c(1,1), distm, NULL, complexity = 0)$C
min(eigen(kdistm)$value)
invkdistm <- solve(kdistm)
plot(x, kdistm[,1], type="l")
plot(x, invkdistm[,1], type="l")
plot(x, abs(invkdistm[,1]), type="l")
plot(x[x>3.8], abs(invkdistm[,1])[x>3.8], type="l")

source("../R/HMC-functions.R")
x <- seq(0,5,0.05)
distm <- as.matrix(dist(x))
kdistm <- calCovCompact1(c(1,1), distm, NULL, complexity = 0)$C
min(eigen(kdistm)$value)
invkdistm <- solve(kdistm)
plot(x, kdistm[,1], type="l")
plot(x, invkdistm[,1], type="l")
plot(x, abs(invkdistm[,1]), type="l")
plot(x[x>3.8], abs(invkdistm[,1])[x>3.8], type="l")

x <- seq(0,10,0.05)
distm <- as.matrix(dist(x))
kdistm <- calCov2(c(1,1), distm, NULL, complexity = 0)$C
min(eigen(kdistm)$value)
invkdistm <- solve(kdistm)
plot(x, kdistm[,1], type="l")
plot(x, invkdistm[,1], type="l")
plot(x, abs(invkdistm[,1]), type="l")
plot(x[x>1], abs(invkdistm[,1])[x>1], type="l")
plot(x[x>2], abs(invkdistm[,1])[x>2], type="l")
plot(x[x>3], abs(invkdistm[,1])[x>3], type="l")
plot(x[x>4], abs(invkdistm[,1])[x>4], type="l")
plot(x[x>5], abs(invkdistm[,1])[x>5], type="l")
plot(x[x>6], abs(invkdistm[,1])[x>6], type="l")

x <- seq(0,6,0.05)
distm <- as.matrix(dist(x))
kdistm <- calCovRBF(c(1,1), distm, NULL, complexity = 0)$C
min(eigen(kdistm)$value)
diag(kdistm) <- diag(kdistm) + abs(min(eigen(kdistm)$value))*2
invkdistm <- solve(kdistm)
plot(x, kdistm[,1], type="l")
plot(x, invkdistm[,1], type="l")
plot(x, abs(invkdistm[,1]), type="l")
plot(x[x>1], abs(invkdistm[,1])[x>1], type="l")
plot(x[x>2], abs(invkdistm[,1])[x>2], type="l")
plot(x[x>3], abs(invkdistm[,1])[x>3], type="l")
plot(x[x>4], abs(invkdistm[,1])[x>4], type="l")
plot(x[x>5], abs(invkdistm[,1])[x>5], type="l")
