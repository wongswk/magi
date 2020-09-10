testthat::context("HMC")
testthat::test_that("HMC runs without error", {
  magi:::hmcTest()
})

testthat::test_that("HMC for normal distribution is correct", {
  out.all.c <- matrix(nrow=1e4, ncol=4)
  out.all.c[1,] <- rep(0,4)
  for(i in 2:nrow(out.all.c)){
    out.normal <- magi:::hmcNormal(out.all.c[i-1,], rep(0.05,4), -Inf, +Inf, 20, TRUE)
    out.all.c[i,] <- out.normal$final
  }
  for(j in 1:4){
    suppressWarnings(checkoutput <- ks.test(out.all.c[,j], "pnorm"))
    testthat::expect_gt(checkoutput$p.value, 1e-5)
  }
})

testthat::test_that("HMC for truncated normal distribution is correct", {
  out.all.c <- matrix(nrow=1e4, ncol=4)
  out.all.c[1,] <- rep(0,4)
  for(i in 2:nrow(out.all.c)){
    out.normal <- magi:::hmcNormal(out.all.c[i-1,], rep(0.05,4), -1, 2, 20, TRUE)
    out.all.c[i,] <- out.normal$final
  }
  ptruncnorm <- function(x){
    out <- (pnorm(x)-pnorm(-1))/(pnorm(2)-pnorm(-1))
    pmin(pmax(0, out), 1)
  }
  for(j in 1:4){
    suppressWarnings(checkoutput <- ks.test(out.all.c[,j], "ptruncnorm"))
    testthat::expect_gt(checkoutput$p.value, 1e-5)
  }
})

testthat::test_that("HMC for generic distribution is correct", {
  llk <- function(x) {
    value = -0.5*sum(x^2)
    gradient = -x
    list(value=value, gradient=gradient)
  }
  hmc_out <- magi::basic_hmcRcpp(llk, c(0.2, 0.2), c(0.1, 0.1), c(-10, -10), c(Inf, Inf), 200, TRUE)
  testthat::expect_equal(llk(hmc_out$final)$value, hmc_out$lprvalue, scale=abs(hmc_out$lprvalue))
})
