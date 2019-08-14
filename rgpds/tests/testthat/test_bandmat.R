testthat::context("band matrix likelihood")
testthat::test_that("band matrix likelihood runs without error", {
  out <- gpds:::bandTest(system.file("testdata/data_band.txt",package="gpds"))
  testthat::expect_equal(out, -655203.16255410481244, tolerance = 1e-3)
})

