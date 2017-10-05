testthat::context("band matrix likelihood")
testthat::test_that("band matrix likelihood runs without error", {
  gpds:::bandTest(system.file("testdata/data_band.txt",package="gpds"))
})

testthat::test_that("band matrix likelihood wrapped runs correctly", {
  datainput <- scan(system.file("testdata/data_band.txt",package="gpds"), 
                    sep = "\n", what = character())
  datainput <- strsplit(datainput, "\t")
  datainput <- lapply(datainput, function(x) as.numeric(na.omit(as.numeric(x))))
  
  foo <- xthetallik(xtheta = datainput[[1]], 
                    Vmphi = datainput[[2]], VKinv = datainput[[3]], VCinv = datainput[[4]],
                    Rmphi = datainput[[5]], RKinv = datainput[[6]], RCinv = datainput[[7]],
                    bandsize = datainput[[8]], n = datainput[[9]], sigma = datainput[[10]],
                    yobs = datainput[[11]])
  outsum <- as.numeric(foo)+sum(attr(foo, "grad"))
  testthat::expect_equal(outsum, -655203.16255410481244)
})

