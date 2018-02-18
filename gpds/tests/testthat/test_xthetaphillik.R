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
dataInputFullObs <- dataInput
dataInputWithMissing <- dataInput
dataInputWithMissing[-seq(1,nrow(dataInputWithMissing),4),] <- NA

xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2])
thetaTest <- pram.true$abc
phiTest <- cbind(pram.true$vphi, pram.true$rphi)
sigmaTest <- rep(pram.true$sigma, 2)

testthat::test_that("xthetaphisigmallik differs to xthetallik and loglikOrig by constant fixing phi sigma", {
  curCovV <- calCov(phiTest[,1], r, signr, kerneltype = "generalMatern")
  curCovR <- calCov(phiTest[,2], r, signr, kerneltype = "generalMatern")
  
  for(dummy in 1:4){
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
      constDiff12 <- -29.84029
      constDiff23 <- -45.94867
    }else{
      dataInput <- dataInputFullObs
      constDiff12 <- -77.772002
      constDiff23 <- -60.65
    }
    
    xthInit <- c(xlatentTest, thetaTest)
    
    out1 <- gpds::xthetallikC(dataInput, curCovV, curCovR, sigmaTest, xthInit)
    out2 <- loglikOrig(xlatentTest,
                       thetaTest,
                       phiTest,
                       sigmaTest[1],
                       dataInput,
                       r,
                       signr,
                       kerneltype = "generalMatern")
    out3 <- xthetaphisigmallikRcpp(xlatentTest,
                                   thetaTest,
                                   phiTest,
                                   sigmaTest,
                                   dataInput,
                                   fn.sim$time)
    
    expect_equal(out1$value - as.numeric(out2), constDiff12, tolerance = 0.00001)
    expect_equal(out3$value - as.numeric(out2), constDiff23, tolerance = 0.001)
  }
})

testthat::test_that("xthetaphisigmallik differs to loglikOrig by constant (the pi part)", {
  for(dummy in 1:4){
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    phiTest <- cbind(pram.true$vphi, pram.true$rphi) * rexp(4)
    sigmaTest <- rep(pram.true$sigma * rexp(1), 2)
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
      constDiff23 <- -45.94867
    }else{
      dataInput <- dataInputFullObs
      constDiff23 <- -60.65
    }
    
    xthInit <- c(xlatentTest, thetaTest)
    
    out2 <- loglikOrig(xlatentTest,
                       thetaTest,
                       phiTest,
                       sigmaTest[1],
                       dataInput,
                       r,
                       signr,
                       kerneltype = "generalMatern")
    out3 <- xthetaphisigmallikRcpp(xlatentTest,
                                   thetaTest,
                                   phiTest,
                                   sigmaTest,
                                   dataInput,
                                   fn.sim$time)
    
    expect_equal(out3$value - as.numeric(out2), constDiff23, tolerance = 0.001)
  }
})


testthat::test_that("xthetaphisigmallik derivatives", {
  
  
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2])
  thetaTest <- pram.true$abc
  phiTest <- cbind(pram.true$vphi, pram.true$rphi)
  sigmaTest <- rep(pram.true$sigma, 2)
  
  
  out <- xthetaphisigmallikRcpp(xlatentTest,
                                thetaTest,
                                phiTest,
                                sigmaTest,
                                dataInputWithMissing,
                                fn.sim$time)
  out$value
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (gpds::xthetallikC(dataInputWithMissing, curCovV, curCovR, cursigma, xthInit1,
                         useBand = FALSE)$value -
         gpds::xthetallikC(dataInputWithMissing, curCovV, curCovR, cursigma, xthInit,
                           useBand = FALSE)$value)/delta
  }
  x <- (gradNum - out$grad)/abs(out$grad)
  testthat::expect_true(all(abs(x) < 1e-3)) # gradient is self-consistent
})


if(interactive()){
  # speed benchmark
  mb <- microbenchmark::microbenchmark(
  gpds::xthetallikC(dataInput, curCovV, curCovR, pram.true$sigma, xthInit),
  loglikOrig(fn.true[seq(1,nrow(fn.true), length=nobs),],
             pram.true$abc,
             c(pram.true$vphi, pram.true$rphi),
             pram.true$sigma,
             dataInput,
             r,
             signr,
             kerneltype = "generalMatern"),
  xthetaphisigmallikRcpp(data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]),
                         pram.true$abc,
                         cbind(pram.true$vphi, pram.true$rphi),
                         rep(pram.true$sigma, 2),
                         dataInput,
                         fn.sim$time)
  )
  print(mb)
}
