library(testthat)
library(magi)

context("x theta sigma log likelihood")

nobs <- 11
noise <- 0.05

VRtrue <- read.csv(system.file("testdata/FN.csv", package="magi"))
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
curCovV$tvecCovInput = fn.sim$time
curCovR <- calCov(phiTest[,2], r, signr, kerneltype = "generalMatern")
curCovR$tvecCovInput = fn.sim$time

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
    
    out2 <- magi:::loglikWithNormalizingConstants(xlatentTest,
                                           thetaTest,
                                           phiTest,
                                           sigmaTest[1],
                                           dataInput,
                                           r,
                                           signr,
                                           kerneltype = "generalMatern")
    out3 <- magi:::xthetasigmallikRcpp(xlatentTest,
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
    
    out2 <- magi:::loglikWithNormalizingConstants(xlatentTest,
                                           thetaTest,
                                           phiTest,
                                           sigmaTest[1],
                                           dataInput,
                                           r,
                                           signr,
                                           kerneltype = "generalMatern")
    out3 <- magi:::xthetasigmallikRcpp(xlatentTest,
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
  priorTemperatureTest <- rexp(3)
  
  out <- magi:::xthetasigmallikRcpp(xlatentTest,
                              thetaTest,
                              sigmaTest,
                              dataInput,
                              list(curCovV, curCovR),
                              priorTemperatureInput=priorTemperatureTest)
  out$value
  
  delta <- 1e-6
  
  # xlatent
  gradNum <- c()
  for(it in 1:length(xlatentTest)){
    xlatentTest1 <- xlatentTest
    xlatentTest1[it] <- xlatentTest1[it] + delta
    gradNum[it] <- 
      (magi:::xthetasigmallikRcpp(xlatentTest1,
                           thetaTest,
                           sigmaTest,
                           dataInput,
                           list(curCovV, curCovR),
                           priorTemperatureInput=priorTemperatureTest)$value -
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
      (magi:::xthetasigmallikRcpp(xlatentTest,
                           thetaTest1,
                           sigmaTest,
                           dataInput,
                           list(curCovV, curCovR),
                           priorTemperatureInput=priorTemperatureTest)$value -
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
      (magi:::xthetasigmallikRcpp(xlatentTest,
                           thetaTest,
                           sigmaTest1,
                           dataInput,
                           list(curCovV, curCovR),
                           priorTemperatureInput=priorTemperatureTest)$value -
         out$value)/delta
  }
  x <- (gradNum - tail(out$grad, 2))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # one combined sigma
  sigmaTest1 <- sigmaTest[1]
  sigmaTest1 <- sigmaTest1 + delta
  out1sigma <- magi:::xthetasigmallikRcpp(xlatentTest,
                                   thetaTest,
                                   sigmaTest[1],
                                   dataInput,
                                   list(curCovV, curCovR),
                                   priorTemperatureInput=priorTemperatureTest)
    
  gradNum <- 
    (magi:::xthetasigmallikRcpp(xlatentTest,
                         thetaTest,
                         sigmaTest1,
                         dataInput,
                         list(curCovV, curCovR),
                         priorTemperatureInput=priorTemperatureTest)$value -
       out1sigma$value)/delta
  
  x <- (gradNum - tail(out1sigma$grad, 1))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
})

testthat::test_that("xthetasigmallik with band approx or mean component", {
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  
  dataInput <- dataInputWithMissing
  constDiff23 <- 34.08103
  
  xthInit <- c(xlatentTest, thetaTest)
  
  outlist <- list()
  for(useMean in c(FALSE, TRUE))
    for(useBand in c(FALSE, TRUE))
      outlist <- c(outlist, list(
        magi:::xthetasigmallikRcpp(xlatentTest,
                            thetaTest,
                            sigmaTest,
                            dataInput,
                            list(curCovV, curCovR),
                            useBand = useBand,
                            useMean = useMean)
      ))
  
  expect_true(all(sapply(outlist, function(x) x$value) == outlist[[1]]$value))
  expect_true(all(sapply(outlist, function(x) x$grad[,1]) == outlist[[1]]$grad[,1]))
})


test_that("xthetasigma sampler can run", {
  stepsize <- rep(0.01, length(c(xlatentTest, thetaTest, sigmaTest)))
  ret <- magi:::xthetasigmaSample(dataInputWithMissing,
                    list(curCovV, curCovR),
                    sigmaTest,
                    c(xlatentTest, thetaTest),
                    stepsize,
                    20,
                    modelName = "FN")
  testthat::expect_equal(length(ret$final), 27)

  ret <- magi:::xthetasigmaSample(dataInputWithMissing,
                    list(curCovV, curCovR),
                    sigmaTest,
                    c(xlatentTest, thetaTest),
                    stepsize,
                    20,
                    modelName = "FN")
  testthat::expect_equal(length(ret$final), 27)
  
  stepsize <- rep(0.001, length(c(xlatentTest, thetaTest, sigmaTest[1])))
  ret <- magi:::xthetasigmaSample(dataInputWithMissing,
                    list(curCovV, curCovR),
                    sigmaTest[1],
                    c(xlatentTest, thetaTest),
                    stepsize,
                    20,
                    loglikflag = "withmeanBand",
                    modelName = "FN")
  testthat::expect_equal(length(ret$final), 26)
})

test_that("xthetasigma sampler can run for Hes1", {
  xsim.obs <- data.frame( 
    time = seq(0, 240, 24),
    X1 = c(0.72, 4.03, 9.55, 6.52, 4.89, 2, 4.27, 7.53, 7.28, 5.25, 1.56), 
    X2 = c(2.3, 2.4, 1.63, 0.97, 0.58, 1.65, 2.33, 1.54, 0.86, 0.26, 1.16), 
    X3 = c(1.4, 3.41, 1.63, 0.11, 5.3, 13, 1.47, 4.66, -0.52, 2.92, 7.73)
  )
  xsim <- setDiscretization(xsim.obs, 1)
  xInit <- data.matrix(xsim[,-1])
  xInit[] <- 0
  thetaInit <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
  sigmaInit <- c(0.8, 0.2, 1.6)
  curphi <- structure(c(24.87, 55.75, 2.52, 60.64, 9.51, 82.61), .Dim = 2:3)
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], abs(outer(xsim$time, xsim$time, '-')), -sign(outer(xsim$time, xsim$time, '-')), 
                      bandsize=20, kerneltype="generalMatern")
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  stepSize <- rep(1e-4, length(xInit) + length(thetaInit) + length(sigmaInit))
  ret <- magi:::xthetasigmaSample(data.matrix(xsim[,-1]), curCov, sigmaInit, c(xInit, thetaInit),
                    stepSize, 20, F, loglikflag = "withmeanBand",
                    priorTemperature = 1, modelName = "Hes1")
  testthat::expect_equal(length(ret$final), 73)
})
