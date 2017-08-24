setwd("~/Workspace/DynamicSys/dynamic-systems/test/")
library(Rcpp)
sourceCpp("../src/wrapper.cpp")

VRtrue <- read.csv("../data/FN.csv")

pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157)
)

nobs <- 6


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
phisigllikTest( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
phisigllik( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), fn.sim[,1:2], TRUE)

phisigSample(data.matrix(fn.sim[,1:2]), r, c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise),
             rep(0.1,5), 20, T)
