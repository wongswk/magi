load("../test/test_xthetaSample.rda")
sourceCpp("../src/wrapper.cpp")
source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

out <- xthetallikTest(dataInput, curCovV, curCovR, cursigma, xthInit)
testthat::expect_equal(out, outExpected)

plot(1/curCovV$Keigen1over)
truncEigen(1/curCovV$Keigen1over)
plot(curCovV$KeigenVec[,1])
plot(curCovV$KeigenVec[,2])
plot(curCovV$KeigenVec[,3])
plot(curCovV$KeigenVec[,4])
plot(curCovV$KeigenVec[,201])

approxCovV <- truncCovByEigen(curCovV, 
                              truncEigen(rev(curCovV$Ceigen1over), 0.85),
                              truncEigen(rev(curCovV$Keigen1over), 0.85))
approxCovV$Keigen1over

x <- approxCovV$KeigenVec%*%diag(approxCovV$Keigen1over)%*%t(approxCovV$KeigenVec)
sum((x - approxCovV$Kinv)^2)/sum(approxCovV$Kinv^2)
diag(abs(x - approxCovV$Kinv)/abs(approxCovV$Kinv))

approxCovR <- truncCovByEigen(curCovR, 
                              truncEigen(rev(curCovR$Ceigen1over), 0.85),
                              truncEigen(rev(curCovR$Keigen1over), 0.85))

outApprox <- xthetallikTest(dataInput, approxCovV, approxCovR, cursigma, xthInit)
outApprox$value - outExpected$value

delta <- 1e-8
gradNum <- c()
for(it in 1:length(xthInit)){
  xthInit1 <- xthInit
  xthInit1[it] <- xthInit1[it] + delta
  gradNum[it] <- 
    (xthetallikTest(dataInput, approxCovV, approxCovR, cursigma, xthInit1)$value -
       xthetallikTest(dataInput, approxCovV, approxCovR, cursigma, xthInit)$value)/delta
}
x <- (gradNum - outApprox$grad)/abs(outApprox$grad)
plot(x)
summary(x)

testthat::expect_equal(outApprox, outExpected)

microbenchmark::microbenchmark(
  xthetallikTest(dataInput, approxCovV, approxCovR, cursigma, xthInit),
  xthetallikTest(dataInput, curCovV, curCovR, cursigma, xthInit)  
)


