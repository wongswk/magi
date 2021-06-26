library(testthat)
context("check pass by reference for R gpcov list")

tvec.nobs <- seq(0,5,0.1)
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
signr <- -sign(foo)
curCovV <- calCov(c(1,1), r, signr)
test_that("curCovV value changed when call c++ test function", {
  magi:::changeGPcovFromC(curCovV)
  
  expect_true(all(curCovV$Cinv==1))
  expect_true(all(curCovV$mphi==2))
  expect_true(all(curCovV$Kinv==3))
  expect_true(all(curCovV$CinvBand==4))
  expect_true(all(curCovV$mphiBand==5))
  expect_true(all(curCovV$KinvBand==6))
  expect_true(all(curCovV$mu==77))
  expect_true(all(curCovV$dotmu==666))
})


test_that("some ad hoc testing on syntax", {
  Amat <- matrix(rnorm(4),2)
  Bmat <- matrix(rnorm(4),2)
  Alist <- list(Cinv=Amat, Brand=Bmat)
  Alist$Cinv <- Alist$Cinv+1

  magi:::cov_r2cpp_t1(Alist)
  magi:::cov_r2cpp_t2(Alist$Cinv)
  expect_equal(Alist$Cinv[1, 1], 0)

  magi:::cov_r2cpp_t2(Amat)
  magi:::cov_r2cpp_t2(Alist$Cinv)
  expect_equal(Amat[1, 1], 0)

  magi:::cov_r2cpp_t3(Bmat)
  expect_equal(Bmat[1, 1], 0)
})
