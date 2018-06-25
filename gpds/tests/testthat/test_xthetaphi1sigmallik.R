library(testthat)
library(gpds)

context("x theta sigma phi1 log likelihood")

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
phi1Test <- phiTest[1,]
sigmaTest <- rep(pram.true$sigma, 2)

curCovV <- calCov(phiTest[,1], r, signr, kerneltype = "generalMatern")
curCovR <- calCov(phiTest[,2], r, signr, kerneltype = "generalMatern")
nophi1Cov <- list(
  calCov(c(1, pram.true$vphi[2]), r, signr, kerneltype = "generalMatern"),
  calCov(c(1, pram.true$rphi[2]), r, signr, kerneltype = "generalMatern")
)

testthat::test_that("xthetaphi1sigmallik reduces to xthetasigmallik", {

  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  
  xthInit <- c(xlatentTest, thetaTest)
  
  out3 <- xthetasigmallikRcpp(xlatentTest,
                              thetaTest,
                              sigmaTest,
                              dataInput,
                              list(curCovV, curCovR))
  out0 <- xthetaphi1sigmallikRcpp(xlatentTest,
                                  thetaTest,
                                  phi1Test,
                                  sigmaTest,
                                  dataInput,
                                  nophi1Cov)
  expect_equal(out3$value / out0$value, 1, tolerance=1e-5)
  expect_equal(out3$grad[1:length(xlatentTest)], 
               out0$grad[1:length(xlatentTest)],
               tolerance=1e-5)
  expect_equal(out3$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))], 
               out0$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))],
               tolerance=1e-5)
  expect_equal(as.numeric(tail(out3$grad, length(sigmaTest))), 
               as.numeric(tail(out0$grad, length(sigmaTest))))
})

testthat::test_that("xthetaphi1sigmallik reduces to xthetasigmallik, even with temperature and mean", {
  priorTemperature = rexp(2)
  nophi1Cov[[1]]$mu <- curCovV$mu <- sin(1:nrow(fn.sim))
  nophi1Cov[[1]]$dotmu <- curCovV$dotmu <- cos(1:nrow(fn.sim))
  nophi1Cov[[2]]$mu <- curCovR$mu <- cos(1:nrow(fn.sim))
  nophi1Cov[[2]]$dotmu <- curCovR$dotmu <- -sin(1:nrow(fn.sim))
  
  for(useMean in c(FALSE, TRUE)){
    for(useBand in c(FALSE, TRUE)){
      xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
      thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
      
      xthInit <- c(xlatentTest, thetaTest)
      
      out3 <- xthetasigmallikRcpp(xlatentTest,
                                  thetaTest,
                                  sigmaTest,
                                  dataInput,
                                  list(curCovV, curCovR),
                                  priorTemperatureInput = priorTemperature,
                                  useBand=useBand,
                                  useMean=useMean)
      out0 <- xthetaphi1sigmallikRcpp(xlatentTest,
                                      thetaTest,
                                      phi1Test,
                                      sigmaTest,
                                      dataInput,
                                      nophi1Cov,
                                      priorTemperatureInput = priorTemperature,
                                      useBand=useBand,
                                      useMean=useMean)
      expect_equal(out3$value / out0$value, 1, tolerance=1e-5)
      expect_equal(out3$grad[1:length(xlatentTest)], 
                   out0$grad[1:length(xlatentTest)],
                   tolerance=1e-5)
      expect_equal(out3$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))], 
                   out0$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))],
                   tolerance=1e-5)
      expect_equal(as.numeric(tail(out3$grad, length(sigmaTest))), 
                   as.numeric(tail(out0$grad, length(sigmaTest))))
    }
  }
})


testthat::test_that("xthetaphi1sigmallik derivatives, withmean useBand", {
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
  phi1Test <- c(pram.true$vphi[1], pram.true$rphi[1]) * exp(rnorm(2))
  
  nophi1Cov[[1]]$mu <- curCovV$mu <- sin(1:nrow(fn.sim))
  nophi1Cov[[1]]$dotmu <- curCovV$dotmu <- cos(1:nrow(fn.sim))
  nophi1Cov[[2]]$mu <- curCovR$mu <- cos(1:nrow(fn.sim))
  nophi1Cov[[2]]$dotmu <- curCovR$dotmu <- -sin(1:nrow(fn.sim))
  
  out <<- xthetaphi1sigmallikRcpp(xlatentTest,
                                  thetaTest,
                                  phi1Test,
                                  sigmaTest,
                                  dataInput,
                                  nophi1Cov,
                                  useMean = TRUE,
                                  useBand = TRUE)
  
  delta <- 1e-6
  
  # xlatent
  gradNum <- c()
  for(it in 1:length(xlatentTest)){
    xlatentTest1 <- xlatentTest
    xlatentTest1[it] <- xlatentTest1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest1,
                               thetaTest,
                               phi1Test,
                               sigmaTest,
                               dataInput,
                               nophi1Cov,
                               useMean = TRUE,
                               useBand = TRUE)$value -
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
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest1,
                               phi1Test,
                               sigmaTest,
                               dataInput,
                               nophi1Cov,
                               useMean = TRUE,
                               useBand = TRUE)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # phi1
  gradNum <- c()
  for(it in 1:length(phi1Test)){
    phi1Test1 <- phi1Test
    phi1Test1[it] <- phi1Test1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest,
                               phi1Test1,
                               sigmaTest,
                               dataInput,
                               nophi1Cov,
                               useMean = TRUE,
                               useBand = TRUE)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+length(thetaTest)+1):(length(xlatentTest)+length(thetaTest)+length(phi1Test))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  
  # two seperate sigma
  gradNum <- c()
  for(it in 1:length(sigmaTest)){
    sigmaTest1 <- sigmaTest
    sigmaTest1[it] <- sigmaTest1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest,
                               phi1Test,
                               sigmaTest1,
                               dataInput,
                               nophi1Cov,
                               useMean = TRUE,
                               useBand = TRUE)$value -
         out$value)/delta
  }
  x <- (gradNum - tail(out$grad, 2))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # one combined sigma
  sigmaTest1 <- sigmaTest[1]
  sigmaTest1 <- sigmaTest1 + delta
  out1sigma <-xthetaphi1sigmallikRcpp(xlatentTest,
                                      thetaTest,
                                      phi1Test,
                                      sigmaTest[1],
                                      dataInput,
                                      nophi1Cov,
                                      useMean = TRUE,
                                      useBand = TRUE)
  
  gradNum <- 
    (xthetaphi1sigmallikRcpp(xlatentTest,
                             thetaTest,
                             phi1Test,
                             sigmaTest1,
                             dataInput,
                             nophi1Cov,
                             useMean = TRUE,
                             useBand = TRUE)$value -
       out1sigma$value)/delta
  
  x <- (gradNum - tail(out1sigma$grad, 1))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
})


testthat::test_that("xthetaphi1sigmallik derivatives", {
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
  phi1Test <- c(pram.true$vphi[1], pram.true$rphi[1]) * exp(rnorm(2))
  
  out <<- xthetaphi1sigmallikRcpp(xlatentTest,
                                  thetaTest,
                                  phi1Test,
                                  sigmaTest,
                                  dataInput,
                                  nophi1Cov)
  
  delta <- 1e-6
  
  # xlatent
  gradNum <- c()
  for(it in 1:length(xlatentTest)){
    xlatentTest1 <- xlatentTest
    xlatentTest1[it] <- xlatentTest1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest1,
                               thetaTest,
                               phi1Test,
                               sigmaTest,
                               dataInput,
                               nophi1Cov)$value -
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
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest1,
                               phi1Test,
                               sigmaTest,
                               dataInput,
                               nophi1Cov)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # phi1
  gradNum <- c()
  for(it in 1:length(phi1Test)){
    phi1Test1 <- phi1Test
    phi1Test1[it] <- phi1Test1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest,
                               phi1Test1,
                               sigmaTest,
                               dataInput,
                               nophi1Cov)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+length(thetaTest)+1):(length(xlatentTest)+length(thetaTest)+length(phi1Test))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  
  # two seperate sigma
  gradNum <- c()
  for(it in 1:length(sigmaTest)){
    sigmaTest1 <- sigmaTest
    sigmaTest1[it] <- sigmaTest1[it] + delta
    gradNum[it] <- 
      (xthetaphi1sigmallikRcpp(xlatentTest,
                               thetaTest,
                               phi1Test,
                               sigmaTest1,
                               dataInput,
                               nophi1Cov)$value -
         out$value)/delta
  }
  x <- (gradNum - tail(out$grad, 2))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # one combined sigma
  sigmaTest1 <- sigmaTest[1]
  sigmaTest1 <- sigmaTest1 + delta
  out1sigma <-xthetaphi1sigmallikRcpp(xlatentTest,
                                      thetaTest,
                                      phi1Test,
                                      sigmaTest[1],
                                      dataInput,
                                      nophi1Cov)
  
  gradNum <- 
    (xthetaphi1sigmallikRcpp(xlatentTest,
                             thetaTest,
                             phi1Test,
                             sigmaTest1,
                             dataInput,
                             nophi1Cov)$value -
       out1sigma$value)/delta
  
  x <- (gradNum - tail(out1sigma$grad, 1))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
})


test_that("xthetaphi1sigmaSample sampler can run", {
  stepsize <- rep(0.01, length(c(xlatentTest, thetaTest, phi1Test, sigmaTest)))
  xthetaphi1sigmaSample(dataInputWithMissing,
                        nophi1Cov,
                        phi1Test,
                        sigmaTest,
                        c(xlatentTest, thetaTest),
                        stepsize,
                        20,
                        modelName = "FN")
  
  xthetaphi1sigmaSample(dataInputWithMissing,
                        nophi1Cov,
                        phi1Test,
                        sigmaTest,
                        c(xlatentTest, thetaTest),
                        stepsize,
                        20,
                        loglikflag = "withmeanBand",
                        modelName = "FN")
  
  stepsize <- rep(0.0001, length(c(xlatentTest, thetaTest, phi1Test, sigmaTest[1])))
  xthetaphi1sigmaSample(dataInputWithMissing,
                        nophi1Cov,
                        phi1Test,
                        sigmaTest[1],
                        c(xlatentTest, thetaTest),
                        stepsize,
                        20,
                        loglikflag = "withmeanBand",
                        modelName = "FN")
})

test_that("xthetaphi1sigmaSample sampler can run for Hes1", {
  xsim.obs <- data.frame(
    time = seq(0, 240, 24),
    X1 = c(0.72, 4.03, 9.55, 6.52, 4.89, 2, 4.27, 7.53, 7.28, 5.25, 1.56),
    X2 = c(2.3, 2.4, 1.63, 0.97, 0.58, 1.65, 2.33, 1.54, 0.86, 0.26, 1.16),
    X3 = c(1.4, 3.41, 1.63, 0.11, 5.3, 13, 1.47, 4.66, -0.52, 2.92, 7.73)
  )
  xsim <- insertNaN(xsim.obs, 1)
  xInit <- data.matrix(xsim[,-1])
  xInit[] <- 0
  thetaInit <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
  sigmaInit <- c(0.8, 0.2, 1.6)
  curphi <- structure(c(24.87, 55.75, 2.52, 60.64, 9.51, 82.61), .Dim = 2:3)
  phi1Init <- curphi[1,]
  curphi[1,] <- 1
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], abs(outer(xsim$time, xsim$time, '-')), -sign(outer(xsim$time, xsim$time, '-')),
                      bandsize=20, kerneltype="generalMatern")
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  stepSize <- rep(1e-4, length(xInit) + length(thetaInit) + length(phi1Init) + length(sigmaInit))
  xthetaphi1sigmaSample(data.matrix(xsim[,-1]), curCov, phi1Init, sigmaInit, c(xInit, thetaInit),
                        stepSize, 20, F, loglikflag = "withmeanBand",
                        priorTemperature = 1, modelName = "Hes1")
})
