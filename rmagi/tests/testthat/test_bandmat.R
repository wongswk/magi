testthat::context("band matrix likelihood")
testthat::test_that("band matrix likelihood runs without error", {
  out <- magi:::bandTest(system.file("testdata/data_band.txt",package="magi"))
  testthat::expect_equal(out, -655203.16255410481244, tolerance = 1e-3)
})

