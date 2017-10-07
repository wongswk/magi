setwd("~/Workspace/DynamicSys/dynamic-systems/test/")
library(Rcpp)
sourceCpp("../src/wrapper.cpp")

VRtrue <- read.csv("../data/FN.csv")

pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157)
)

nobs <- 41


library(parallel)
source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

#noise level
noise <- 0.5
pram.true$sigma <- noise

fn.true <- VRtrue
fn.true$time <- seq(0,20,0.05)
fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]

tvec.nobs <- fn.sim$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

n.iter <- 100  # number of HMC iterations
phisig <- matrix(NA,n.iter,5)   # phi and sigma

phisig[1,] <- rep(1,5)

##### Reference values (truth)
bestCovV <- calCov( c( 1.9840824, 1.1185157) )
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
logliknoODE.mar( bestCovV, bestCovR, noise, fn.sim[,1:2])
xc <- phisigllikTest( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
xr <- phisigllik( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), fn.sim[,1:2], TRUE)

x0 <- c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise)
phisigllik(x0, fn.sim[,1:2], TRUE)
for(i in 1:5){
  # check derivative
  x1 = x0
  x1[i] = x1[i] + 1e-9
  print(as.numeric((phisigllik(x1, fn.sim[,1:2]) - phisigllik(x0, fn.sim[,1:2]))/1e-9))
}

phisigllik( c(-1, -1, -1, -1, -1), fn.sim[,1:2], TRUE)
phisigllikTest( c(-1, -1, -1, -1, -1), data.matrix(fn.sim[,1:2]), r)
# need to handle contraint on HMC

testpoint <- abs(rnorm(5))
xc <- phisigllikTest( testpoint, data.matrix(fn.sim[,1:2]), r)
xr <- phisigllik( testpoint, fn.sim[,1:2], TRUE)
xc$value - xr
xc$grad - attr(xr, "grad")

phisigSample(data.matrix(fn.sim[,1:2]), r, c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise),
             rep(0.03,5), 20, F)


xthetaSample(data.matrix(fn.sim[,1:2]), bestCovV, bestCovR, noise, 
             rep(1, nobs*2+3), rep(0.03,nobs*2+3), 1, T)

xthetallik(matrix(1, nrow=nobs, ncol=2), rep(1,3), bestCovV, bestCovR, noise, data.matrix(fn.sim[,1:2]), grad = T)

#### test RBF kernel ####
testpoint <- abs(rnorm(5))
xc <- phisigllikTest( testpoint, data.matrix(fn.sim[,1:2]), r, "rbf")
xr <- phisigllik( testpoint, fn.sim[,1:2], TRUE, "rbf")
as.numeric(xc$value - xr)
xc$grad - attr(xr, "grad")

x0 <- c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise)
phisigllikTest(x0, data.matrix(fn.sim[,1:2]), r, "rbf")$grad
for(i in 1:5){
  # check derivative
  x1 = x0
  x1[i] = x1[i] + 1e-9
  print(as.numeric((phisigllikTest(x1, data.matrix(fn.sim[,1:2]), r, "rbf")$value - 
                      phisigllikTest(x0, data.matrix(fn.sim[,1:2]), r, "rbf")$value)/1e-9))
}

phisigSample(data.matrix(fn.sim[,1:2]), r, c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise),
             rep(0.03,5), 20, F, "rbf")

#### test compact1 kernel ####
testpoint <- abs(rnorm(5))
xc <- phisigllikTest( testpoint, data.matrix(fn.sim[,1:2]), r, "compact1")
xr <- phisigllik( testpoint, fn.sim[,1:2], TRUE, "compact1")
as.numeric(xc$value - xr)
xc$grad - attr(xr, "grad")

x0 <- c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise)
phisigllikTest(x0, data.matrix(fn.sim[,1:2]), r, "rbf")$grad
for(i in 1:5){
  # check derivative
  x1 = x0
  x1[i] = x1[i] + 1e-9
  print(as.numeric((phisigllikTest(x1, data.matrix(fn.sim[,1:2]), r, "rbf")$value - 
                      phisigllikTest(x0, data.matrix(fn.sim[,1:2]), r, "rbf")$value)/1e-9))
}

phisigSample(data.matrix(fn.sim[,1:2]), r, c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise),
             rep(0.03,5), 20, F, "rbf")
