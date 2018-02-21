library(testthat)
library(gpds)

context("x theta sigma log likelihood")

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

curCovV <- calCov(phiTest[,1], r, signr, kerneltype = "generalMatern")
curCovR <- calCov(phiTest[,2], r, signr, kerneltype = "generalMatern")

testthat::test_that("xthetasigmallik differs to xthetallik and loglikOrig by constant fixing phi sigma", {
  
  realDiff <- sapply(1:40, function(dummy){ 
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
      constDiff23 <- 34.08103
    }else{
      dataInput <- dataInputFullObs
      constDiff23 <- 48.78405
    }
    
    xthInit <- c(xlatentTest, thetaTest)
    
    out2 <- loglikWithNormalizingConstants(xlatentTest,
                                           thetaTest,
                                           phiTest,
                                           sigmaTest[1],
                                           dataInput,
                                           r,
                                           signr,
                                           kerneltype = "generalMatern")
    out3 <- xthetasigmallikRcpp(xlatentTest,
                                thetaTest,
                                sigmaTest,
                                dataInput,
                                list(curCovV, curCovR))
    
    out3$value - as.numeric(out2) - constDiff23
  })
  expect_gt(mean(realDiff < 1e-3), 0.95)
})

testthat::test_that("xthetasigmallik differs to loglikOrig by constant (the pi part)", {
  # phi could give numerical instability issue
  realDiff <- sapply(1:40, function(dummy){
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
      constDiff23 <- 34.08103
    }else{
      dataInput <- dataInputFullObs
      constDiff23 <- 48.78405
    }
    
    xthInit <- c(xlatentTest, thetaTest)
    
    out2 <- loglikWithNormalizingConstants(xlatentTest,
                                           thetaTest,
                                           phiTest,
                                           sigmaTest[1],
                                           dataInput,
                                           r,
                                           signr,
                                           kerneltype = "generalMatern")
    out3 <- xthetasigmallikRcpp(xlatentTest,
                                thetaTest,
                                sigmaTest,
                                dataInput,
                                list(curCovV, curCovR))
    out3$value - as.numeric(out2) - constDiff23
  })
  expect_gt(mean(abs(realDiff) < 0.001), 0.4)
})


testthat::test_that("xthetasigmallik derivatives", {
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
  
  out <<- xthetasigmallikRcpp(xlatentTest,
                              thetaTest,
                              sigmaTest,
                              dataInput,
                              list(curCovV, curCovR))
  out$value
  
  delta <- 1e-6
  
  # xlatent
  gradNum <- c()
  for(it in 1:length(xlatentTest)){
    xlatentTest1 <- xlatentTest
    xlatentTest1[it] <- xlatentTest1[it] + delta
    gradNum[it] <- 
      (xthetasigmallikRcpp(xlatentTest1,
                           thetaTest,
                           sigmaTest,
                           dataInput,
                           list(curCovV, curCovR))$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[1:length(xlatentTest)])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # theta
  gradNum <- c()
  for(it in 1:length(thetaTest)){
    thetaTest1 <- thetaTest
    thetaTest1[it] <- thetaTest1[it] + delta
    gradNum[it] <- 
      (xthetasigmallikRcpp(xlatentTest,
                           thetaTest1,
                           sigmaTest,
                           dataInput,
                           list(curCovV, curCovR))$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  
  # two seperate sigma
  gradNum <- c()
  for(it in 1:length(sigmaTest)){
    sigmaTest1 <- sigmaTest
    sigmaTest1[it] <- sigmaTest1[it] + delta
    gradNum[it] <- 
      (xthetasigmallikRcpp(xlatentTest,
                           thetaTest,
                           sigmaTest1,
                           dataInput,
                           list(curCovV, curCovR))$value -
         out$value)/delta
  }
  x <- (gradNum - tail(out$grad, 2))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # one combined sigma
  sigmaTest1 <- sigmaTest[1]
  sigmaTest1 <- sigmaTest1 + delta
  out1sigma <- xthetasigmallikRcpp(xlatentTest,
                                   thetaTest,
                                   sigmaTest[1],
                                   dataInput,
                                   list(curCovV, curCovR))
    
  gradNum <- 
    (xthetasigmallikRcpp(xlatentTest,
                         thetaTest,
                         sigmaTest1,
                         dataInput,
                         list(curCovV, curCovR))$value -
       out1sigma$value)/delta
  
  x <- (gradNum - tail(out1sigma$grad, 1))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
})

test_that("xthetasigma sampler can run", {
  stepsize <- rep(0.01, length(c(xlatentTest, thetaTest, sigmaTest)))
  xthetasigmaSample(dataInputWithMissing,
                    list(curCovV, curCovR),
                    sigmaTest,
                    c(xlatentTest, thetaTest),
                    stepsize,
                    20,
                    modelName = "FN")
  
  stepsize <- rep(0.001, length(c(xlatentTest, thetaTest, sigmaTest[1])))
  xthetasigmaSample(dataInputWithMissing,
                    list(curCovV, curCovR),
                    sigmaTest[1],
                    c(xlatentTest, thetaTest),
                    stepsize,
                    20,
                    modelName = "FN")
})
