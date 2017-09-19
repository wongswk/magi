load("../diag/replicate-HMC-leave-mode.rda")
sourceCpp("../src/wrapper.cpp")

# xth.formal[1,] has llik around -70 using the rescaled llik.
#[1] -70.44782
#attr(,"components")
#          [,1]       [,2]      [,3]
#[1,] -14.26526 -10.663875 -12.07067
#[2,] -15.88596  -8.829259  -8.73280

burnin <- 100
stepLow.scaler <- c()
accepts <- c()
accepts[1] <- 0
lliklist <- c()
lliklist[1] <- -70.44782
n.iter <- 100
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, 20, T)
  xth.formal[t,] <- foo$final
  accepts[t] <- foo$acc
  stepLow.scaler[t] <- mean(stepLow)
  
  if (t < burnin) {
    if (mean(tail(accepts[1:t],100)) > 0.9) {
      stepLow <- stepLow * 1.01
    } else if (mean(tail(accepts[1:t],100)) < 0.6) {
      stepLow <- stepLow * .99
    }
  }
  lliklist[t] <- foo$lpr
  
  if( t %% 10 == 0) show(c(t, mean(tail(accepts[1:t],100)), foo$lpr, tail(foo$final,3)))
}


