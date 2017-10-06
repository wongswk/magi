library(testthat)

testthat::test_that("on attach messages", {
  expect_message(library(gpds))
})

test_check("gpds")

#' add test for other kernels
#' add test for drawing from GP smoothing
#' check test code coverage