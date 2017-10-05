testthat::context("HMC")
testthat::test_that("HMC runs without error", {
  gpds:::hmcTest()
})

testthat::test_that("HMC for normal distribution is correct", {
  out.all.c <- matrix(nrow=1e4, ncol=4)
  out.all.c[1,] <- rep(0,4)
  for(i in 2:nrow(out.all.c)){
    out.normal <- gpds:::hmcNormal(out.all.c[i-1,], rep(0.05,4), -Inf, +Inf, 20, TRUE)
    out.all.c[i,] <- out.normal$final
  }
  for(j in 1:4){
    testthat::expect_gt(ks.test(out.all.c[,4], "pnorm")$p.value, 1e-5)
  }
})

testthat::test_that("HMC for truncated normal distribution is correct", {
  out.all.c <- matrix(nrow=1e4, ncol=4)
  out.all.c[1,] <- rep(0,4)
  for(i in 2:nrow(out.all.c)){
    out.normal <- gpds:::hmcNormal(out.all.c[i-1,], rep(0.05,4), -1, 2, 20, TRUE)
    out.all.c[i,] <- out.normal$final
  }
  ptruncnorm <- function(x){
    out <- (pnorm(x)-pnorm(-1))/(pnorm(2)-pnorm(-1))
    pmin(pmax(0, out), 1)
  }
  for(j in 1:4){
    testthat::expect_gt(ks.test(out.all.c[,4], "ptruncnorm")$p.value, 1e-5)
  }
})
