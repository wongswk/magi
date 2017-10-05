testthat::context("parallel tempering")
testthat::test_that("parallel tempering runs without error", {
  gpds:::paralleltemperingTest1()
  gpds:::paralleltemperingTest2()
})