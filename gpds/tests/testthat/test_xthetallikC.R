testthat::context("xthetallikC")

load(system.file("testdata/test_xthetaSample.rda", package="gpds"))

testthat::test_that("xthetallikC runs without error and is correct", {
  out <- gpds::xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit)
  testthat::expect_equal(out, outExpected)
})

testthat::test_that("examine low rank approximation", {
  plot(1/curCovV$Keigen1over)
  truncEigen(1/curCovV$Keigen1over)
  plot(curCovV$KeigenVec[,1])
  plot(curCovV$KeigenVec[,2])
  plot(curCovV$KeigenVec[,3])
  plot(curCovV$KeigenVec[,4])
  plot(curCovV$KeigenVec[,201])
  
  approxCovVgood <- truncCovByEigen(curCovV, 
                                    truncEigen(rev(curCovV$Ceigen1over), 0.99),
                                    truncEigen(rev(curCovV$Keigen1over), 0.99))
  approxCovRgood <- truncCovByEigen(curCovR, 
                                    truncEigen(rev(curCovR$Ceigen1over), 0.99),
                                    truncEigen(rev(curCovR$Keigen1over), 0.99))
  
  outApprox <- xthetallikC(dataInput, 
                           approxCovVgood, 
                           approxCovRgood, 
                           cursigma, xthInit)
  testthat::expect_lt((outApprox$value-outExpected$value)/outExpected$value, 1e-2)
  
  approxCovV <- truncCovByEigen(curCovV, 
                                truncEigen(rev(curCovV$Ceigen1over), 0.85),
                                truncEigen(rev(curCovV$Keigen1over), 0.85))
  
  x <- approxCovV$KeigenVec%*%diag(approxCovV$Keigen1over)%*%t(approxCovV$KeigenVec)
  sum((x - approxCovV$Kinv)^2)/sum(approxCovV$Kinv^2)
  diag(abs(x - approxCovV$Kinv)/abs(approxCovV$Kinv))
  
  approxCovR <- truncCovByEigen(curCovR, 
                                truncEigen(rev(curCovR$Ceigen1over), 0.85),
                                truncEigen(rev(curCovR$Keigen1over), 0.85))
  
  outApprox <- xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit)
  abs((outApprox$value - outExpected$value)/outExpected$value)
  
  delta <- 1e-8
  gradNum <- c()
  for(it in 1:length(xthInit)){
    xthInit1 <- xthInit
    xthInit1[it] <- xthInit1[it] + delta
    gradNum[it] <- 
      (xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit1)$value -
         xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit)$value)/delta
  }
  x <- (gradNum - outApprox$grad)/abs(outApprox$grad)
  testthat::expect_true(all(abs(x) < 1e-3))
  
  x <- microbenchmark::microbenchmark(
    xthetallikC(dataInput, approxCovV, approxCovR, cursigma, xthInit),
    xthetallikC(dataInput, approxCovVgood, approxCovRgood, cursigma, xthInit),
    xthetallikC(dataInput, curCovV, curCovR, cursigma, xthInit),
    times=10
  )
  x <- tapply(x$time, x$expr, mean)
  testthat::expect_true(x[1] < x[2] & x[2] < x[3])
})



