library(testthat)

testthat::test_that("on attach messages", {
  expect_message(library(gpds))
})

test_check("gpds")

#' check test code coverage