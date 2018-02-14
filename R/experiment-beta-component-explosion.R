library(testthat)
library(gpds)

config <- list(
  nobs = 26,
  noise = 0.1,
  kernel = "generalMatern",
  seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
  npostplot = 50,
  loglikflag = "withmeanBand",
  filllevel = 5
)
explosion_result <- list()

for(filllevel in 0:5){
  config$filllevel <- filllevel
  
  pram.true <- list(
    theta=c(0.2,0.2,3),
    x0=c(-1, 1),
    rphi=c(0.9486433, 3.2682434),
    vphi=c(1.9840824, 1.1185157),
    sigma=config$noise
  )
  
  modelODE <- function(t, state, parameters) {
    list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
  }
  times <- seq(0,20,0.1)
  xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
  xtrue <- data.frame(xtrue)
  if(filllevel==0){
    matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)  
  }
  
  xtrueFunc <- lapply(2:ncol(xtrue), function(j)
    approxfun(xtrue[, "time"], xtrue[, j]))
  
  xsim <- data.frame(xtrue)
  
  set.seed(config$seed)
  for(j in 1:(ncol(xsim)-1)){
    xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise)  
  }
  
  xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
  
  xsim <- insertNaN(xsim.obs,config$filllevel)
  
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
  
  
  fn <- function(par) -phisigllikC( par, data.matrix(xsim.obs[,-1]), 
                                    r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(xsim.obs[,-1]), 
                                              r.nobs, config$kernel)$grad)
  marlikmap <- optim(rep(1, 5), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma <- marlikmap$par[5]
  curphi[] <- marlikmap$par[1:4]
  
  cursigma
  curphi
  
  curCov <- lapply(1:(ncol(xsim.obs)-1), function(j){
    covEach <- calCov(curphi[, j], r, signr, bandsize=config$bandsize, 
                      kerneltype=config$kernel)
    covEach$mu[] <- mean(xsim.obs[,j+1])
    covEach
  })
  
  
  
  xInit <- sapply(xtrueFunc, function(f) f(xsim$time))
  
  xthetallikOurC <- xthetallikC( data.matrix(xsim[,-1]),
                                 curCov[[1]],
                                 curCov[[2]],
                                 cursigma,
                                 c(xInit, pram.true$theta))
  
  
  xthetallikOurR <- xthetallik( xInit,
                                pram.true$theta,
                                curCov[[1]],
                                curCov[[2]],
                                cursigma,
                                data.matrix(xsim[,-1]),
                                T)
  
  testthat::expect_equal(as.numeric(xthetallikOurR), xthetallikOurC$value, tolerance = 1e-5)
  testthat::expect_equal(attr(xthetallikOurR, "grad"), as.numeric(xthetallikOurC$grad), 
                         tolerance = 1e-5)
  
  if(filllevel==0){
    curphi0 <- curphi
    cursigma0 <- cursigma
    print(curphi)
    print(cursigma)
  }else{
    expect_equal(curphi0, curphi)
    expect_equal(cursigma0, cursigma)
  }
  explosion_result[[as.character(filllevel)]] <- list(
    nrow(xsim),
    attr(xthetallikOurR, "components")
  )
}

explosion_result_component <- sapply(explosion_result, function(x) x[[2]], simplify = "array")
dimnames(explosion_result_component)[[1]] <- c("V","R")
dimnames(explosion_result_component)[[2]] <- c("y-x","f-mx K f-mx","x C x")
matplot(t(log(-explosion_result_component[,"x C x",])), typ="b",
        xlab="number of discretization", ylab="log quadratic form for x", xaxt='n')
axis(1, 1:6, sapply(explosion_result, function(x) x[[1]]))
title(expression(x^T * C^-1 * x))
mtext(paste("random seed =", config$seed))
matplot(t(log(-explosion_result_component[,"f-mx K f-mx",])), type="b",
        xlab="number of discretization", ylab="log quadratic form for dot-x", xaxt='n')
axis(1, 1:6, sapply(explosion_result, function(x) x[[1]]))
title(expression((f - m*x)^T * K^-1 * (f - m*x)))
mtext(paste("random seed =", config$seed))
