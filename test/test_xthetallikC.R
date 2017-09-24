load("../test/test_xthetaSample.rda")
sourceCpp("../src/wrapper.cpp")
out <- xthetallikTest(dataInput, curCovV, curCovR, cursigma, xthInit)
testthat::expect_equal(out, outExpected)