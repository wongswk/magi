library(magi)

config <- list(
  nobs = 11,
  noise = c(4,1,8)*0.2,
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  bandsize = 20,
  hmcSteps = 500,
  n.iter = 5000,
  burninRatio = 0.50,
  stepSizeFactor = 1,
  filllevel = 3,
  modelName = "hes1"
)


config$ndis <- (config$nobs-1)*2^config$filllevel+1

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(0.5, 2, 1),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::hes1modelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
#matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- xtrue

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- setDiscretization(xsim.obs,config$filllevel)

tvec.full <- xsim$time
tvec.nobs <- xsim.obs$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -magi:::phisigllikC( par, data.matrix(xsim.obs[,1+j]),
                                    r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(magi:::phisigllikC( par, data.matrix(xsim.obs[,1+j]),
                                              r.nobs, config$kernel)$grad)
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}
cursigma
curphi

j <- 1
# plot(xsim.obs[, "time"], data.matrix(xsim.obs[,1+j]))
magi:::phisigllikC(c(curphi[,j], cursigma[j]), data.matrix(xsim.obs[,1+j]),
            r.nobs, config$kernel)
magi:::phisigllikC(c(100, 100, cursigma[j]), data.matrix(xsim.obs[,1+j]),
            r.nobs, config$kernel)


xtrueAtDiscretization <- sapply(xtrueFunc, function(f) f(xsim[,"time"]))
test_that("xthetaphisigmallikRcpp can run", {
  out <- xthetaphisigmallikRcpp( xtrueAtDiscretization,
                          pram.true$theta,
                          matrix(pram.true$phi, nrow = 2),
                          pram.true$sigma,
                          data.matrix(xsim[,-1]),
                          xsim$time,
                          "Hes1")
  testthat::expect_equal(length(out$grad), 259)
})
