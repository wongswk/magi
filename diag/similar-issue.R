#### tempered mixture distribution example ####
#' original sample from a mixture where the peak dominates
#' then heat up system (rescale likelihood), we can see sampler leaves the peak
#' and loglik goes down as more samples are drawn
source("../R/helper/basic_hmc.R")
source("../R/helper/utilities.r")

mixden <- function(x, grad=F){
  # log(0.05*dnorm(x,-1)+0.9*dnorm(x,0)+0.05*dnorm(x,1))
  sharpsd <- 0.01
  temperature <- 0.1
  pp <- 0.01
  ret <- log(pp*dnorm(x,0,sharpsd)+(1-pp)*dnorm(x,0,1))/temperature
  if(grad){
    val <- (pp*dnorm(x,0,sharpsd)+(1-pp)*dnorm(x,0,1))
    gradient <- pp*dnorm(x,0,sharpsd)*(-x/sharpsd^2) + (1-pp)*dnorm(x,0,1)*(-x/1^2)
    gradient <- gradient/val/temperature
    attr(ret, "gradient") <- gradient
  }
  ret
}

mixden.rescale <- function(x, grad=F){
  ret <- mixden(x, grad=grad)
  ret <- ret*26/201
  if(grad){
    attr(ret, "gradient") <- attr(ret, "gradient")*26/201
  }
  ret
}

mixden(0.001,T)
(mixden(0.001+1e-9)-mixden(0.001))/1e-9

mixden.rescale(0.001,T)
(mixden.rescale(0.001+1e-9)-mixden.rescale(0.001))/1e-9

mixden.rescale(0.1,T)
(mixden.rescale(0.1+1e-9)-mixden.rescale(0.1))/1e-9


plot.function(function(x) exp(mixden(x)), from=-0.2,to=0.2,n=1e4)
plot.function(function(x) exp(mixden.rescale(x)), from=-3,to=3,n=1e4)

n.iter <- 1e4
xsam.orig <- 0
for(t in 2:n.iter){
  foo <- basic_hmc(mixden, step=runif(1,0.01,0.02), nsteps= 20, initial=xsam.orig[t-1], return.traj = T)
  xsam.orig[t] <- foo$final
}
plot(xsam.orig)
plot(mixden.rescale(xsam.orig))
hist(xsam.orig, breaks=20)

xsam.rescale <- 0
for(t in 2:n.iter){
  foo <- basic_hmc(mixden.rescale, step=runif(1,0.1,0.2), nsteps= 20, initial=xsam.rescale[t-1], return.traj = T)
  xsam.rescale[t] <- foo$final
}
plot(mixden.rescale(xsam.rescale))
hist(xsam.rescale, breaks=200)

