library(testthat)
library(gpds)

context("x theta phi log likelihood")

nobs <- 11
noise <- 0.05

VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=noise
)
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


dataInput <- data.matrix(fn.sim[,1:2])

curCovV <- calCov(pram.true$vphi, r, signr, kerneltype = "generalMatern")
curCovR <- calCov(pram.true$rphi, r, signr, kerneltype = "generalMatern")

fn.true[,1:2] <- fn.true[,1:2] * rexp(length(fn.true[,1:2]))
pram.true$abc <- pram.true$abc * rexp(length(pram.true$abc))
xthInit <- c(data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]), pram.true$abc)
out1 <- gpds::xthetallikC(dataInput, curCovV, curCovR, pram.true$sigma, xthInit)
out2 <- loglikOrig(fn.true[seq(1,nrow(fn.true), length=nobs),],
                   pram.true$abc,
                   c(pram.true$vphi, pram.true$rphi),
                   pram.true$sigma,
                   dataInput,
                   r,
                   signr,
                   kerneltype = "generalMatern")

expect_equal(out1$value - as.numeric(out2), -77.772002)

out3 <- xthetaphisigmallikRcpp(data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]),
                               pram.true$abc,
                               cbind(pram.true$vphi, pram.true$rphi),
                               pram.true$sigma,
                               dataInput,
                               fn.sim$time)
