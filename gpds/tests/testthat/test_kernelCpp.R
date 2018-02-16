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

