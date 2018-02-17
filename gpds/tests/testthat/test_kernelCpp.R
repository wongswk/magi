library(testthat)
testthat::context("testing GP kernel cpp functions")

phitrue <- list(
  generalMatern = c(2.04, 1.313)
)

xtime <- seq(0,20,0.5)

kerneltype <- "generalMatern"
testthat::test_that("check C, dCdphi for generalMatern",{
  
  egcovR <- calCov(phitrue[[kerneltype]], 
                   as.matrix(dist(xtime)),
                   -sign(outer(xtime,xtime,'-')),
                   kerneltype = kerneltype)
  dimnames(egcovR$Cdoubleprime) <- NULL
  dimnames(egcovR$C) <- NULL
  
  egcovC <- generalMaternCovRcpp(phitrue[[kerneltype]], outer(xtime,xtime,'-'))
  expect_equal(egcovC$C, egcovR$C, check.attributes = FALSE)
  expect_equal(egcovC$dCdphiCube, egcovR$dCdphiCube, check.attributes = FALSE)
  expect_equal(egcovC$Cprime, egcovR$Cprime, check.attributes = FALSE)
  expect_equal(egcovC$Cdoubleprime, egcovR$Cdoubleprime, check.attributes = FALSE)
  
  if(interactive()){
    mb <- microbenchmark::microbenchmark(
      egcovR <- calCovGeneralMatern(phitrue[[kerneltype]], 
                                    as.matrix(dist(xtime)),
                                    -sign(outer(xtime,xtime,'-'))),
      egcovC <- generalMaternCovRcpp(phitrue[[kerneltype]], outer(xtime,xtime,'-'))
    )
    print(mb)
  }
  
})

testthat::test_that("underflow for generalMatern when used in hes1 with loocvllik",{
  xsim.obs <- data.frame(time = c(0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240), 
                         X1 = c(-0.661, 5.074, 8.82, 7.028, 4.869, 2.527, 4.498, 8.262, 9.443, 4.876, 1.561), 
                         X2 = c(2.5, 2.443, 1.365, 0.941, 0.827, 1.565, 2.442, 1.417, 1.099, 0.599, 1.589), 
                         X3 = c(0.162, 1.802, 1.725, -0.101, 2.804, 17.063, 5.93, -1.825, 0.816, 1.362, 10.697))
  foo <- outer(xsim.obs$time, t(xsim.obs$time),'-')[,1,]
  r.nobs <- abs(foo)
  
  out <- phisigloocvllikC(  c(117.3426722,   0.1,   0.9607398), data.matrix(xsim.obs$X3), 
                            r.nobs, "generalMatern")
  
  expect_true(is.finite(out$value))
  expect_true(all(is.finite(out$grad)))
  
  egcovC <- generalMaternCovRcpp(c(117.3426722, 0.1), foo)
  
})
