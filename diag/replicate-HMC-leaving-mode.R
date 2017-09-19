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
lliklist[1] <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                            xth.formal[1,], rstep*100, 1, T, F)$lpr
n.iter <- 100
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, 20, T, F)
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

layout(matrix(1:6,ncol=2,byrow = T))
hist(xth.formal[1:n.iter,405])
plot(xth.formal[1:n.iter,405], main="c")
hist(xth.formal[1:n.iter,404])
plot(xth.formal[1:n.iter,404], main="b")
hist(xth.formal[1:n.iter,403])
plot(xth.formal[1:n.iter,403], main="a")


lliklist_rescaled <- sapply(1:n.iter, function(t)
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xth.formal[t,], rstep*100, 1, T, T)$lpr)
layout(1:2)
plot(lliklist, main="lliklist (unscaled)")
plot(lliklist_rescaled, main="lliklist (rescaled)")

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
lliklist[1] <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                            xth.formal[1,], rstep*100, 1, T, T)$lpr
n.iter <- 100
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, 20, T, T)
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
layout(matrix(1:6,ncol=2,byrow = T))
hist(xth.formal[1:n.iter,405])
plot(xth.formal[1:n.iter,405], main="c")
hist(xth.formal[1:n.iter,404])
plot(xth.formal[1:n.iter,404], main="b")
hist(xth.formal[1:n.iter,403])
plot(xth.formal[1:n.iter,403], main="a")

lliklist_unscaled <- sapply(1:n.iter, function(t)
  xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
               xth.formal[t,], rstep*100, 1, T, F)$lpr)

layout(1:2)
plot(lliklist, main="lliklist (rescaled)")
plot(lliklist_unscaled, main="lliklist (unscaled)")
