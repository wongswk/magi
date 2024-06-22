library(magi)

test_that("output custom object gpcov", {
    out = magi:::maternCovTestOutput(1:2, cbind(1:2, 2:1))
    expect_equal(out$C, cbind(c(0.8286492424, 0.523994108), c(0.5239941088, 0.82864924)))
})

input_cov = magi:::maternCovTestOutput(1:2, cbind(1:2, 2:1))
input_cov$Cinv = solve(input_cov$C)
input_cov$mphi = input_cov$C + 1
input_cov$Kinv = input_cov$C + 2
input_cov$CinvBand = input_cov$C + 3
input_cov$mphiBand = input_cov$C + 4
input_cov$KinvBand = input_cov$C + 5
input_cov$mu = c(11, 12)
input_cov$dotmu = c(21, 22)
input_cov$bandsize = 2

test_that("input custom object gpcov", {
    skip("workaround Rcpp issue #1308")
    expect_equal(magi:::maternCovTestInput(input_cov), input_cov$Cinv)
})
