testthat::context("band matrix likelihood")
testthat::test_that("band matrix likelihood runs without error", {
  gpds:::bandTest(system.file("testdata/data_band.txt",package="gpds"))
})