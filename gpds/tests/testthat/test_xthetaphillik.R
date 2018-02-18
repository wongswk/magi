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
 
  realDiff <- sapply(1:40, function(dummy){ 
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
      constDiff12 <- 16.10664
    }else{
      dataInput <- dataInputFullObs
      constDiff12 <- -17.12206
    }
    
    xthInit <- c(xlatentTest, thetaTest)
    
    out1 <- gpds::xthetallikC(dataInput, curCovV, curCovR, sigmaTest, xthInit)
    out2 <- loglikWithNormalizingConstants(xlatentTest,
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
    
    c(out1$value - as.numeric(out2) - constDiff12, 
      out3$value - as.numeric(out2) )
  })
  expect_gt(mean(realDiff < 1e-3), 0.95)
})

testthat::test_that("xthetaphisigmallik differs to loglikOrig by constant (the pi part)", {
  # phi could give numerical instability issue
  realDiff <- sapply(1:40, function(dummy){
    xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
    thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
    phiTest <- cbind(pram.true$vphi, pram.true$rphi) * exp(rnorm(4))
    sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
    
    if(dummy %% 2 == 0){
      dataInput <- dataInputWithMissing
    }else{
      dataInput <- dataInputFullObs
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
    out3 <- xthetaphisigmallikRcpp(xlatentTest,
                                   thetaTest,
                                   phiTest,
                                   sigmaTest,
                                   dataInput,
                                   fn.sim$time)
    out3$value - as.numeric(out2)
  })
  expect_gt(mean(abs(realDiff) < 0.01), 0.4)
})


testthat::test_that("xthetaphisigmallik derivatives", {
  xlatentTest <- data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]) * rexp(length(fn.true[,1:2]))
  thetaTest <- pram.true$abc * rexp(length(pram.true$abc))
  phiTest <- cbind(pram.true$vphi, pram.true$rphi) * exp(rnorm(4))
  sigmaTest <- rep(pram.true$sigma * exp(rnorm(1)), 2)
  
  out <- xthetaphisigmallikRcpp(xlatentTest,
                                thetaTest,
                                phiTest,
                                sigmaTest,
                                dataInputWithMissing,
                                fn.sim$time)
  out$value
  
  delta <- 1e-5
  
  # xlatent
  gradNum <- c()
  for(it in 1:length(xlatentTest)){
    xlatentTest1 <- xlatentTest
    xlatentTest1[it] <- xlatentTest1[it] + delta
    gradNum[it] <- 
      (xthetaphisigmallikRcpp(xlatentTest1,
                              thetaTest,
                              phiTest,
                              sigmaTest,
                              dataInputWithMissing,
                              fn.sim$time)$value -
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
      (xthetaphisigmallikRcpp(xlatentTest,
                              thetaTest1,
                              phiTest,
                              sigmaTest,
                              dataInputWithMissing,
                              fn.sim$time)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+1):(length(xlatentTest)+length(thetaTest))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # phi
  gradNum <- c()
  for(it in 1:length(phiTest)){
    phiTest1 <- phiTest
    phiTest1[it] <- phiTest1[it] + delta
    gradNum[it] <- 
      (xthetaphisigmallikRcpp(xlatentTest,
                              thetaTest,
                              phiTest1,
                              sigmaTest,
                              dataInputWithMissing,
                              fn.sim$time)$value -
         out$value)/delta
  }
  x <- (gradNum - out$grad[(length(xlatentTest)+length(thetaTest)+1):
                             (length(xlatentTest)+length(thetaTest)+length(phiTest))])/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # two seperate sigma
  gradNum <- c()
  for(it in 1:length(sigmaTest)){
    sigmaTest1 <- sigmaTest
    sigmaTest1[it] <- sigmaTest1[it] + delta
    gradNum[it] <- 
      (xthetaphisigmallikRcpp(xlatentTest,
                              thetaTest,
                              phiTest,
                              sigmaTest1,
                              dataInputWithMissing,
                              fn.sim$time)$value -
         out$value)/delta
  }
  x <- (gradNum - tail(out$grad, 2))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
  
  # one combined sigma
  sigmaTest1 <- sigmaTest[1]
  sigmaTest1 <- sigmaTest1 + delta
  out1sigma <- xthetaphisigmallikRcpp(xlatentTest,
                                      thetaTest,
                                      phiTest,
                                      sigmaTest[1],
                                      dataInputWithMissing,
                                      fn.sim$time)
  gradNum <- 
    (xthetaphisigmallikRcpp(xlatentTest,
                            thetaTest,
                            phiTest,
                            sigmaTest1,
                            dataInputWithMissing,
                            fn.sim$time)$value -
       out1sigma$value)/delta
  
  x <- (gradNum - tail(out1sigma$grad, 1))/abs(gradNum)
  testthat::expect_true(all(abs(x) < 5e-3)) # gradient is self-consistent
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
