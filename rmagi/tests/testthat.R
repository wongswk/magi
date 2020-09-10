library(testthat)

testthat::test_that("on attach messages", {
  expect_message(library(magi))
})

test_check("magi")

#' check test code coverage