testthat::context("testing GP kernel functions")

#--- on sampled data, fitted derivative and true derivative are similar ----
#' the variance of derivative seems wrong

xtime <- seq(0,20,0.1)
phitrue <- list(
  compact1 = c(2.618, 6.381, 0.152, 9.636),
  rbf = c(0.838, 0.307, 0.202, 0.653),
  matern = c(2.04, 1.313, 0.793, 3.101),
  periodicMatern = c(2.04, 1.313, 9, 0.793, 3.101, 9),
  generalMatern = c(2.04, 1.313, 0.793, 3.101),
  rationalQuadratic = c(85.278, 5.997, 35.098, 14.908)
)

uprange <- 3.5
downrange <- 0.4

set.seed(123)
for(kerneltype in c("compact1","rbf","matern","periodicMatern","generalMatern",
                    "rationalQuadratic")){
  if(kerneltype=="rationalQuadratic") downrange <- 0.07
  testthat::test_that("check Cprime and Cdoubleprime", {
    xtime <- seq(0,2,0.01)
    delta <- mean(diff(xtime))
    egcov2 <- calCov(head(phitrue[[kerneltype]], length(phitrue[[kerneltype]])/2), 
                     as.matrix(dist(xtime)),
                     -sign(outer(xtime,xtime,'-')),
                     kerneltype = kerneltype)
    
    testthat::expect_equal( (egcov2$Cprime[-1,1] + egcov2$Cprime[-nrow(egcov2$Cprime),1])/2, 
                           diff(egcov2$C[,1])/delta, 
                           tolerance = 0.0001, scale=max(abs(egcov2$Cprime[-1,1])))
    testthat::expect_equal( (egcov2$Cdoubleprime[-1,1] + egcov2$Cdoubleprime[-nrow(egcov2$Cdoubleprime),1])/2, 
                           -diff(egcov2$Cprime[,1])/0.01, 
                           tolerance = 0.001, scale=max(abs(egcov2$Cdoubleprime[-1,1])))
  })
  
  testthat::test_that("check dCdphi", {
    if(kerneltype == "periodicMatern") skip("periodicMatern dCdphi is wrong for the 3rd parameter")
    xtime <- seq(0,2,0.1)
    testpintPhi <- head(phitrue[[kerneltype]], length(phitrue[[kerneltype]])/2)
    delta <- 1e-4
    egcov0 <- calCov(testpintPhi, 
                     as.matrix(dist(xtime)),
                     -sign(outer(xtime,xtime,'-')),
                     kerneltype = kerneltype)
    
    for(it in 1:length(testpintPhi)){
      testpintPhi1 <- testpintPhi
      testpintPhi1[it] <- testpintPhi1[it] + delta
      egcov1 <- calCov(testpintPhi1, 
                       as.matrix(dist(xtime)),
                       -sign(outer(xtime,xtime,'-')),
                       kerneltype = kerneltype)
      expect_equal((egcov1$C - egcov0$C)/delta, egcov0$dCdphiCube[,,it],
                   tolerance = 1e-3, scale = max(abs(egcov0$dCdphiCube[,,it])))
    }
  })
  
  curCovV <- calCov(head(phitrue[[kerneltype]], length(phitrue[[kerneltype]])/2), 
                    as.matrix(dist(xtime)),
                    -sign(outer(xtime,xtime,'-')),
                    kerneltype = kerneltype)
  
  curCovR <- calCov(tail(phitrue[[kerneltype]], length(phitrue[[kerneltype]])/2), 
                    as.matrix(dist(xtime)),
                    -sign(outer(xtime,xtime,'-')),
                    kerneltype = kerneltype)
  
  
  
  fn.true <- read.csv(system.file("testdata/FN.csv", package="magi"))
  fn.true$time <- fn.true$time <- seq(0,20,0.05)
  
  abc <- c(0.2, 0.2, 3)
  
  fn.true$dVtrue = with(fn.true, abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
  fn.true$dRtrue = with(fn.true, -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))
  
  # V
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), curCovV$C)
  plot(fn.true$time, fn.true$Vtrue, lwd=2, type="l")
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8, add=T)
  mtext(paste("V -", kerneltype))
  
  testthat::test_that(paste("level shoule be in the range -",kerneltype), {
    rangeRatio <- range(fn.true$Vtrue) / range(draws)
    expect_lt(max(rangeRatio), uprange)  
    expect_gt(min(rangeRatio), downrange)
  })
  
  x <- draws[1,]
  dxNum <- c(diff(head(x,2))/(xtime[2]-xtime[1]),
             (x[-(1:2)]-x[1:(length(x)-2)])/(xtime[3]-xtime[1]),
             diff(tail(x,2))/(xtime[2]-xtime[1]))
  dxMean <- curCovV$mphi%*%x
  plot(xtime, dxNum, type="l", main="derivative")
  lines(xtime, dxMean, col=2)
  mtext(paste("V -", kerneltype))
  
  test_that(paste("check if derivative variance Kmat too small: V of FN - ", kerneltype),{
    if(kerneltype == "rbf") skip("rbf derivative variance Kmat too small")
    if(kerneltype == "rationalQuadratic") skip("rationalQuadratic derivative variance Kmat too small")
    expect_lt(mean(abs((dxNum-dxMean)/sqrt(diag(curCovR$Kphi)))), 10)
  })
  
  
  testthat::test_that(paste("derivative difference should be small on the eye -",kerneltype), {
    expect_lt(mean(abs((dxNum - dxMean)[2:(length(dxMean)-1)]/diff(range(dxMean)))),0.05)  
  })
  
  # R
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), curCovR$C)
  plot(fn.true$time, fn.true$Rtrue, lwd=2, type="l")
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8, add=T)
  mtext(paste("R -", kerneltype))
  
  testthat::test_that(paste("level shoule be in the range -",kerneltype), {
    rangeRatio <- range(fn.true$Vtrue) / range(draws)
    expect_lt(max(rangeRatio), uprange)  
    expect_gt(min(rangeRatio), downrange)
  })
  
  x <- draws[1,]
  dxNum <- c(diff(head(x,2))/(xtime[2]-xtime[1]),
             (x[-(1:2)]-x[1:(length(x)-2)])/(xtime[3]-xtime[1]),
             diff(tail(x,2))/(xtime[2]-xtime[1]))
  dxMean <- curCovR$mphi%*%x
  plot(xtime, dxNum, type="l", main="derivative")
  lines(xtime, dxMean, col=2)
  mtext(paste("R -", kerneltype))
  
  test_that(paste("check if derivative variance Kmat too small: R of FN - ", kerneltype),{
    if(kerneltype == "rbf") skip("rbf derivative variance Kmat too small")
    expect_lt(mean(abs((dxNum-dxMean)/sqrt(diag(curCovR$Kphi)))), 10)
  })
  
  testthat::test_that(paste("derivative difference should be small on the eye -",kerneltype), {
    expect_lt(mean(abs((dxNum - dxMean)[2:(length(dxMean)-1)]/diff(range(dxMean)))),0.05)  
  })
  
  #--- on real data, fitted derivative and true derivative are quite different ----
  fn.true <- fn.true[seq(1,401,2),]
  
  # V
  x <- fn.true$Vtrue
  Ctrue <- curCovV$C
  dxNum <- fn.true$dVtrue
  dxMean <- curCovV$mphi%*%x
  plot(fn.true$time, dxNum, type="l")
  lines(fn.true$time, dxMean, col=2)
  plot((dxNum - dxMean)/sqrt(diag(curCovV$Kphi)), type="l")
  
  testthat::test_that(paste("derivative difference should be small on the eye - realdata ",kerneltype), {
    expect_lt(mean(abs((dxNum - dxMean)[2:(length(dxMean)-1)]/diff(range(dxMean)))),0.01)  
  })
  
  test_that(paste("check if derivative variance Kmat too small: realdata V of FN - ", kerneltype),{
    expect_lt(mean(abs((dxNum-dxMean)/sqrt(diag(curCovV$Kphi)))), 5)
  })
  
  # R
  x <- fn.true$Rtrue
  dxNum <- fn.true$dRtrue
  dxMean <- curCovR$mphi%*%x
  plot((dxNum - dxMean)/sqrt(diag(curCovR$Kphi)), type="l")
  
  testthat::test_that(paste("derivative difference should be small on the eye - realdata ",kerneltype), {
    expect_lt(mean(abs((dxNum - dxMean)[2:(length(dxMean)-1)]/diff(range(dxMean)))),0.01)  
  })
  
  test_that(paste("check if derivative variance Kmat too small: realdata R of FN - ", kerneltype),{
    if(kerneltype == "rationalQuadratic") skip("rationalQuadratic derivative variance Kmat too small")
    testthat::expect_lt(mean(abs((dxNum-dxMean)/sqrt(diag(curCovR$Kphi)))), 5)
  })
  
  
}

testthat::test_that("linear kernel gives linear curve", {
  xtime <- seq(0,5,0.1)
  egcov <- magi:::calCovLinear(c(1,1), xtime, complexity = 0)
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), egcov$C)
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8)
  mtext("linear kernel")
  diffdraws <- t(diff(t(draws)))
  testthat::expect_true(all(abs(diffdraws - rowMeans(diffdraws)) < 1e-5))
  # can be reparametrized to be translation invariant if SigmaMat is not diagonal  
})

testthat::test_that("Neural Network kernel gives rapid changes around 0", {
  xtime <- seq(-4,4,0.1)
  egcov <- magi:::calCovNeuralNetwork(c(1,1), xtime, complexity = 0)
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), egcov$C)
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8)
  title("Neural Network kernel")
  mtext("not suitable for us: don't know fast moving part a priori")
})

testthat::test_that("Warping Matern kernel with sin/cos", {
  xtime <- seq(-2,5,0.1)
  periodicity <- 2
  
  r <- as.matrix(dist(xtime))
  signr <- -sign(outer(xtime, xtime, "-"))
  egcov2 <- magi:::calCovPeriodicWarpMatern(c(1,1, periodicity), r, signr, complexity = 3)
  
  draws <- MASS::mvrnorm(2, rep(0, length(xtime)), egcov2$C)
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8)
  title("Warping Matern kernel with sin/cos")
})

testthat::test_that("modulated squared exponential kernel", {
  xtime <- seq(-6,6,0.1)
  egcov <- magi:::calCovModulatedRBF(c(1,1,1), xtime, complexity = 0)
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), egcov$C)
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8)
  title("modulated squared exponential kernel")
  mtext("not suitable for us: don't have tightening at both ends")
})

testthat::test_that("Rational Quadratic kernel", {
  xtime <- seq(-2,5,0.1)
  
  r <- as.matrix(dist(xtime))
  signr <- -sign(outer(xtime, xtime, "-"))
  egcov2 <- magi:::calCovRationalQuadratic(c(1,1), r, signr, complexity = 3)
  
  draws <- MASS::mvrnorm(7, rep(0, length(xtime)), egcov2$C)
  matplot(xtime, t(draws), type="l", lty = 2:8, col= 2:8)
  title("Rational Quadratic kernel")
  mtext("could be promising")
})
