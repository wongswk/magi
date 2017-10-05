testthat::context("HMC-functions")
library(testthat)
library(gpds)

#### set up: taken from HMC-v3 ####
nobs <- 41
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

set.seed(123)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]

tvec.nobs <- fn.sim$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

testthat::test_that("xthetallikC runs without error and is correct", {
  phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
  fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[,1:2]), r)$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[,1:2]), r)$grad)
  marlikmap <<- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
  testthat::expect_equal(marlikmap$par,
                         c(1.94957912838359, 1.06073179913222, 0.703951749929781, 2.74496621777115, 
                           0.0497280144187114))
})

testthat::test_that("calCov runs without error and is correct", {
  curCovV <<- calCov(marlikmap$par[1:2], r, signr)
  curCovR <<- calCov(marlikmap$par[3:4], r, signr)
  varnames <- c("C", "Cprime", "Cdoubleprime", "Cinv", "mphi", "Kphi", "Kinv")
  curCovV.checksum <- sapply(curCovV[varnames], sum)
  expect_equal(curCovV.checksum, 
               structure(c(387.258476932602, -4.33680868994202e-19, 18.9915012185013, 
                           4.54580185807132, 2.77122075287295e-16, 3.56097233411689, 565.445767848231
               ), .Names = c("C", "Cprime", "Cdoubleprime", "Cinv", "mphi", 
                             "Kphi", "Kinv")))
  curCovR.checksum <- sapply(curCovR[varnames], sum)
  expect_equal(curCovR.checksum, 
               structure(c(335.611085502683, -3.25260651745651e-18, 5.6676295529988, 
                           5.81607236583291, -1.591760143832e-12, 0.0154034834126055, 168699.947542838
               ), .Names = c("C", "Cprime", "Cdoubleprime", "Cinv", "mphi", 
                             "Kphi", "Kinv")))
})
cursigma <- marlikmap$par[5]

testthat::test_that("getMeanCurve runs without error and is correct", {
  startVR <<- rbind(getMeanCurve(fn.sim$time, fn.sim$Vtrue, fn.sim$time, 
                                t(marlikmap$par[1:2]), sigma.mat=matrix(cursigma)),
                   getMeanCurve(fn.sim$time, fn.sim$Rtrue, fn.sim$time, 
                                t(marlikmap$par[3:4]), sigma.mat=matrix(cursigma)))
  testthat::expect_equal(sum(startVR), 16.6934654159305)
})
startVR <- t(startVR)

testthat::test_that("loglikOrig and loglik runs without error and is correct", {
  out <- loglikOrig(fn.true[seq(1,nrow(fn.true), length=nobs),],
                    pram.true$abc,
                    c(pram.true$vphi, pram.true$rphi),
                    pram.true$sigma,
                    fn.sim,
                    r,
                    signr)
  expect_equal(out, 
               structure(454.809052651206, 
                         components = structure(c(104.147795635414, 
                                                  108.143924612949, 18.4183735575236, 
                                                  140.020977657482, -1.5461761617613, 
                                                  85.6241573495981), .Dim = 2:3)))
})

testthat::context("xthetallikC")

dataInput <- data.matrix(fn.sim[,1:2])
xthInit <- c(data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]), pram.true$abc)

outExpectedvalue <- -94.8205825207303

testthat::test_that("xthetallikC runs without error and is correct", {
  out <- gpds::xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit)
  
  testthat::expect_equal(out$value, -94.8205825207303)
  testthat::expect_equal(sum(out$grad), 167.746373733369)
})

testthat::test_that("examine low rank approximation", {
  plot(1/curCovV$Keigen1over)
  truncEigen(1/curCovV$Keigen1over)
  plot(curCovV$KeigenVec[,1])
  plot(curCovV$KeigenVec[,2])
  plot(curCovV$KeigenVec[,3])
  plot(curCovV$KeigenVec[,4])
  
  approxCovVgood <- truncCovByEigen(curCovV, 
                                    truncEigen(rev(curCovV$Ceigen1over), 0.99),
                                    truncEigen(rev(curCovV$Keigen1over), 0.99))
  approxCovRgood <- truncCovByEigen(curCovR, 
                                    truncEigen(rev(curCovR$Ceigen1over), 0.99),
                                    truncEigen(rev(curCovR$Keigen1over), 0.99))
  
  outApprox <- xthetallikC(dataInput, 
                           approxCovVgood, 
                           approxCovRgood, 
                           cursigma, xthInit)
  testthat::expect_lt((outApprox$value-outExpectedvalue)/outExpectedvalue, 1e-2)
  
  approxCovV <- truncCovByEigen(curCovV, 
                                truncEigen(rev(curCovV$Ceigen1over), 0.85),
                                truncEigen(rev(curCovV$Keigen1over), 0.85))
  
  x <- approxCovV$KeigenVec%*%diag(approxCovV$Keigen1over)%*%t(approxCovV$KeigenVec)
  sum((x - approxCovV$Kinv)^2)/sum(approxCovV$Kinv^2)
  diag(abs(x - approxCovV$Kinv)/abs(approxCovV$Kinv))
  
  approxCovR <- truncCovByEigen(curCovR, 
                                truncEigen(rev(curCovR$Ceigen1over), 0.85),
                                truncEigen(rev(curCovR$Keigen1over), 0.85))
  
  outApprox <- xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit)
  abs((outApprox$value - outExpectedvalue)/outExpectedvalue)
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit1)$value -
         xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit)$value)/delta
  }
  x <- (gradNum - outApprox$grad)/abs(outApprox$grad)
  testthat::expect_true(all(abs(x) < 1e-3))
  
  x <- microbenchmark::microbenchmark(
    xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit),
    xthetallikC(dataInput, approxCovVgood, approxCovRgood, cursigma, xthInit),
    xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit),
    times=10
  )
  x <- tapply(x$time, x$expr, mean)
  testthat::expect_true(x[1] < x[3] & x[2] < x[3])
})

bandsize <- 10
testthat::test_that("examine band matrix approximation", {
  curCovVband <<- bandCov(curCovV, bandsize)
  curCovRband <<- bandCov(curCovR, bandsize)
})
