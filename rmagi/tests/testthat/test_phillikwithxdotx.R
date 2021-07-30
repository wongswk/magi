library(testthat)
library(magi)

context("phillikwithxdotx")

times <- seq(0, 60*4, by = 2)
pram.true <- list(
  phi = c(15, 20),
  sigma = 1,
  kernel = "matern"
)

tvec.nobs <- times
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

covSim <- calCov(pram.true$phi, r.nobs, signr.nobs, kerneltype = pram.true$kernel)

x <- MASS::mvrnorm(n=4, mu=rep(0, length(times)), covSim$C)
matplot(times, t(x), type="l")
x <- x[1,]
plot(times, x, type="l")
plot(times, covSim$mphi%*%x, type="l")
dx <- (x[-1] - x[-length(x)])/mean(diff(times))
dx <- apply(cbind(c(NA,dx), c(dx,NA)), 1, mean, na.rm=TRUE)
lines(times, dx, col=2)
plot(dx - covSim$mphi%*%x)

test_that("MLE for phillikwithxdotx", {
  fn <- function(par) -magi:::phillikwithxdotx( par, x, dx,
                                         r.nobs, signr.nobs, pram.true$kernel)
  marlikmap <- optim(rep(1, 2), fn, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf), hessian = TRUE)
  
  expect_gt(
    magi:::phillikwithxdotx( marlikmap$par, x, dx,
                      r.nobs, signr.nobs, pram.true$kernel),
    magi:::phillikwithxdotx( c(pram.true$phi), x, dx,
                      r.nobs, signr.nobs, pram.true$kernel)
  )
  
  expect_lt(
    magi:::phillikwithxdotx( marlikmap$par, x, dx,
                      r.nobs, signr.nobs, pram.true$kernel) - 
      magi:::phillikwithxdotx( c(pram.true$phi), x, dx,
                        r.nobs, signr.nobs, pram.true$kernel),
    1e3
  )
  
  expect_equal(
    mvtnorm::dmvnorm(t(x), sigma = covSim$C, log = TRUE) +
      mvtnorm::dmvnorm(t(dx), mean = covSim$mphi%*%x, sigma = covSim$Kphi, log = TRUE),
    magi:::phillikwithxdotx( c(pram.true$phi), x, dx,
                      r.nobs, signr.nobs, pram.true$kernel),
    tolerance = 1
  )
  
  expect_lt(
    magi:::phillikwithxdotx( marlikmap$par, x, covSim$mphi%*%x,
                      r.nobs, signr.nobs, pram.true$kernel),
    magi:::phillikwithxdotx( c(pram.true$phi), x, covSim$mphi%*%x,
                      r.nobs, signr.nobs, pram.true$kernel)
  )
  
})
