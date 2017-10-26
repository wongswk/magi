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

testthat::test_that("phisigllikC runs without error and is correct", {
  phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
  fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[,1:2]), r)$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[,1:2]), r)$grad)
  marlikmap <<- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
  testthat::expect_equal(marlikmap$par,
                         c(1.94957912838359, 1.06073179913222, 0.703951749929781, 2.74496621777115, 
                           0.0497280144187114),
                         tolerance = 1e-5)
  
  fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[,1:2]), r, "compact1")$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[,1:2]), r, "compact1")$grad)
  marlikmapCompact1 <<- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
  testthat::expect_equal(marlikmapCompact1$par,
                         c(2.04871398302633, 3.59648132314111, 0.625313733996474, 8.96656240950113, 
                           0.0431289806093459),
                         tolerance = 1e-5)
  
  fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[,1:2]), r, "rbf")$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[,1:2]), r, "rbf")$grad)
  marlikmapRbf <<- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
  testthat::expect_equal(marlikmapRbf$par,
                         c(1.58950017432859, 0.594486181493354, 0.452955008218075, 1.57490985649202, 
                           0.0554703243636422),
                         tolerance = 1e-5)
})


testthat::test_that("calCov runs without error and is correct", {
  curCovV <<- calCov(marlikmap$par[1:2], r, signr)
  curCovR <<- calCov(marlikmap$par[3:4], r, signr)
  varnames <- c("C", "Cprime", "Cdoubleprime", "Cinv", "mphi", "Kphi", "Kinv")
  curCovV.checksum <- sapply(curCovV[varnames], function(x) sum(abs(x)))
  outExpect <- c(387.258476932602, 289.658301140539, 369.687853445841, 
                           1859.64654809337, 224.085486080059, 35.0849556047908, 565.445767848267)
  expect_equal(curCovV.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
  curCovR.checksum <- sapply(curCovR[varnames], function(x) sum(abs(x)))
  outExpect <- c(335.611085502683, 96.4763661864788, 45.1452675162526, 
                 443546.019796851, 244.910143264427, 0.154503453415642, 168699.947543871)
  expect_equal(curCovR.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
  
  curCovVcompact1 <<- calCov(marlikmap$par[1:2], r, signr, kerneltype = "compact1")
  curCovRcompact1 <<- calCov(marlikmap$par[3:4], r, signr, kerneltype = "compact1")
  curCovV.checksum <- sapply(curCovVcompact1[varnames], function(x) sum(abs(x)))
  outExpect <- c(115.084361649859, 205.278332612955, 2131.40178197008, 
                 37.6922719587419, 144.786133517796, 2075.23766151204, 2.71838147604252)
  expect_equal(curCovV.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
  curCovR.checksum <- sapply(curCovRcompact1[varnames], function(x) sum(abs(x)))
  outExpect <- c(102.781796009223, 103.966939672376, 187.641691435036, 
                 1951.74861141087, 179.878679502857, 42.9071355730163, 120.568332217125)
  expect_equal(curCovR.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
  
  curCovVrbf <<- calCov(marlikmap$par[1:2], r, signr, kerneltype = "rbf")
  curCovRrbf <<- calCov(marlikmap$par[3:4], r, signr, kerneltype = "rbf")
  curCovV.checksum <- sapply(curCovVrbf[varnames], function(x) sum(abs(x)))
  outExpect <- c(407.840065042279, 293.008616795784, 338.429857457948, 
                 272835661.009244, 736.757284051975, 0.00565329191451157, 374288748.912434)
  expect_equal(curCovV.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
  curCovR.checksum <- sapply(curCovRrbf[varnames], function(x) sum(abs(x)))
  outExpect <- c(354.860844103142, 95.7538793879521, 43.1253728981135, 
                 786566415.192592, 152.053562220004, 5.00881002256369e-05, 708528973.947655)
  expect_equal(curCovR.checksum/outExpect, rep(1, length(outExpect)),
               tolerance = 1e-4, check.attributes = FALSE)
})
cursigma <- marlikmap$par[5]

testthat::test_that("getMeanCurve runs without error and is correct", {
  startVR <<- rbind(getMeanCurve(fn.sim$time, fn.sim$Vtrue, fn.sim$time, 
                                t(marlikmap$par[1:2]), sigma.mat=matrix(cursigma)),
                   getMeanCurve(fn.sim$time, fn.sim$Rtrue, fn.sim$time, 
                                t(marlikmap$par[3:4]), sigma.mat=matrix(cursigma)))
  testthat::expect_equal(sum(startVR), 16.6934654159305, tolerance = 1e-5)
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
                                                  85.6241573495981), .Dim = 2:3)),
               tolerance = 1e-5)
})

testthat::context("xthetallikC")

dataInput <- data.matrix(fn.sim[,1:2])
xthInit <- c(data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]), pram.true$abc)

#### compact1 kernel ####
outExpectedvalue <- -55.73911

testthat::test_that("compact1 - xthetallikC runs without error and is correct", {
  out <- gpds::xthetallikC(dataInput, curCovVcompact1, curCovRcompact1, cursigma, xthInit)
  
  testthat::expect_equal(out$value, outExpectedvalue, tolerance = 1e-5, scale = outExpectedvalue)
  gradExpect <- 143.910609069213
  testthat::expect_equal(sum(out$grad), gradExpect, tolerance = 1e-4, scale = gradExpect)
})

bandsize <- 15
testthat::test_that("compact1 - examine band matrix approximation", {
  curCovVband <<- bandCov(curCovVcompact1, bandsize)
  curCovRband <<- bandCov(curCovRcompact1, bandsize)
  outBandApprox <<- xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)
  testthat::expect_lt(abs((outBandApprox$value - outExpectedvalue)/outExpectedvalue), 1e-3)
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit1)$value -
         xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)$value)/delta
  }
  x <- (gradNum - outBandApprox$grad)/abs(outBandApprox$grad)
  testthat::expect_true(all(abs(x) < 1e-3))
})

#### rbf kernel #### doesn't fit well
outExpectedvalue <- -6017094

testthat::test_that("rbf - xthetallikC runs without error and is correct", {
  out <- gpds::xthetallikC(dataInput, curCovVrbf, curCovRrbf, cursigma, xthInit)
  
  testthat::expect_equal(out$value, outExpectedvalue, tolerance = 1e-5, scale = outExpectedvalue)
  testthat::expect_equal(sum(out$grad), 2781668, tolerance = 1e-1)
})

bandsize <- 15
testthat::test_that("rbf - examine band matrix approximation", {
  curCovVband <<- bandCov(curCovVrbf, bandsize)
  curCovRband <<- bandCov(curCovRrbf, bandsize)
  outBandApprox <<- xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)
  testthat::expect_gt(abs((outBandApprox$value - outExpectedvalue)/outExpectedvalue), 0.25)
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit1)$value -
         xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)$value)/delta
  }
  x <- (gradNum - outBandApprox$grad)/abs(outBandApprox$grad)
  testthat::expect_true(all(abs(x) < 1e-3)) # gradient is self-consistent
})

#### matern kernel ####
outExpectedvalue <- -94.8205825207303

testthat::test_that("xthetallikC runs without error and is correct", {
  out <- gpds::xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit)
  
  testthat::expect_equal(out$value, outExpectedvalue, tolerance = 1e-5, scale = outExpectedvalue)
  gradExpect <- 167.746373733369
  testthat::expect_equal(sum(out$grad), gradExpect, tolerance = 1e-4, scale = gradExpect)
})

testthat::test_that("xthetallik_rescaledC runs without error and compare to non-scaled", {
  out <- gpds::xthetallik_rescaledC(dataInput, curCovV, curCovR, cursigma, xthInit)
  
  testthat::expect_equal(out$value, outExpectedvalue, tolerance = 1e-5, scale = outExpectedvalue)
  gradExpect <- 167.746373733369
  testthat::expect_equal(sum(out$grad), gradExpect, tolerance = 1e-4, scale = gradExpect)
})

testthat::test_that("xthetallik_withmuC runs without error and compare to zero-mean", {
  out <- gpds::xthetallik_withmuC(dataInput, curCovV, curCovR, cursigma, xthInit)
  
  testthat::expect_equal(out$value, outExpectedvalue, tolerance = 1e-5, scale = outExpectedvalue)
  gradExpect <- 167.746373733369
  testthat::expect_equal(sum(out$grad), gradExpect, tolerance = 1e-4, scale = gradExpect)
  
  dotmu <- fODE(pram.true$abc, data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2]))
  
  curCovV_withmu <- curCovV
  curCovR_withmu <- curCovR
  curCovV_withmu$mu <- fn.true[seq(1,nrow(fn.true), length=nobs),1]
  curCovR_withmu$mu <- fn.true[seq(1,nrow(fn.true), length=nobs),2]
  curCovV_withmu$dotmu <- dotmu[,1]
  curCovR_withmu$dotmu <- dotmu[,2]
  
  out <- gpds::xthetallik_withmuC(dataInput, curCovV_withmu, curCovR_withmu, cursigma, xthInit)
  
  dataWN <- dataInput - data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs),1:2])
  xthWN <- c(rep(0,length(xthInit)-3), 0,0,1)
  outWN <- gpds::xthetallikC(dataWN, curCovV, curCovR, cursigma, xthWN)
  
  testthat::expect_equal(out, outWN, tolerance = 1e-5)
})

bandsize <- 15
testthat::test_that("examine band matrix approximation", {
  curCovVband <<- bandCov(curCovV, bandsize)
  curCovRband <<- bandCov(curCovR, bandsize)
  outBandApprox <<- xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)
  testthat::expect_lt(abs((outBandApprox$value - outExpectedvalue)/outExpectedvalue), 1e-3)
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit1)$value -
         xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit)$value)/delta
  }
  x <- (gradNum - outBandApprox$grad)/abs(outBandApprox$grad)
  testthat::expect_true(all(abs(x) < 1e-3))
})


testthat::test_that("examine low rank approximation", {
  plot(1/curCovV$Keigen1over)
  truncEigen(1/curCovV$Keigen1over)
  plot(curCovV$KeigenVec[,1])
  plot(curCovV$KeigenVec[,2])
  plot(curCovV$KeigenVec[,3])
  plot(curCovV$KeigenVec[,4])
  
  approxCovVgood <- truncCovByEigen(curCovV, 
                                    truncEigen(rev(curCovV$Ceigen1over), 0.95),
                                    truncEigen(rev(curCovV$Keigen1over), 0.95))
  approxCovRgood <- truncCovByEigen(curCovR, 
                                    truncEigen(rev(curCovR$Ceigen1over), 0.95),
                                    truncEigen(rev(curCovR$Keigen1over), 0.95))
  
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
    xthetallikBandApproxC(dataInput, curCovVband, curCovRband, cursigma, xthInit),
    xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit),
    times=50
  )
  # print(x)
  x <- tapply(x$time, x$expr, median)
  testthat::expect_true(x[1] < x[3] & x[2] < x[3])
})



testthat::test_that("band matrix likelihood wrapped runs correctly", {
  datainput <- scan(system.file("testdata/data_band.txt",package="gpds"), 
                    sep = "\n", what = character(), quiet=TRUE)
  datainput <- strsplit(datainput, "\t")
  datainput <- lapply(datainput, function(x) as.numeric(na.omit(as.numeric(x))))
  
  covVpart <- curCovVband
  covVpart$mphiBand <- datainput[[2]]
  covVpart$KinvBand <- datainput[[3]]
  covVpart$CinvBand <- datainput[[4]]
  covVpart$bandsize <- datainput[[8]]
  
  covRpart <- curCovRband
  covRpart$mphiBand <- datainput[[5]]
  covRpart$KinvBand <- datainput[[6]]
  covRpart$CinvBand <- datainput[[7]]
  covRpart$bandsize <- datainput[[8]]
  yobs <- matrix(datainput[[11]], ncol=2)
  yobs[yobs==-99999] <- NaN
  foo <- xthetallikBandApproxC(yobs, covVpart, covRpart, datainput[[10]], datainput[[1]])
  
  outsum <- as.numeric(foo$value)+sum(foo$grad)
  testthat::expect_equal(outsum, -655203.16255410481244, tolerance = 1e-3)
})

