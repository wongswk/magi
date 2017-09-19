#### tempered mixture distribution example ####
#' original sample from a mixture where the peak dominates
#' then heat up system (rescale likelihood), we can see sampler leaves the peak
#' and loglik goes down as more samples are drawn
mixden <- function(x){
  # log(0.05*dnorm(x,-1)+0.9*dnorm(x,0)+0.05*dnorm(x,1))
  pp <- 0.00001
  log(pp*dnorm(x,0,0.0001)+(1-pp)*dnorm(x,0,1))/0.1
}

mixden.rescale <- function(x){
  mixden(x)*26/201
}

plot.function(function(x) exp(mixden(x)), from=-0.2,to=0.2,n=1e4)
plot.function(function(x) exp(mixden.rescale(x)), from=-3,to=3,n=1e4)
n.iter <- 1e4
xsam.orig <- 0
for(t in 2:n.iter){
  prop <- rnorm(1,sd=0.0001)+xsam.orig[t-1]
  if(log(runif(1)) < mixden(prop)-mixden(xsam.orig[t-1])){
    xsam.orig[t] <- prop
  }else{
    xsam.orig[t] <- xsam.orig[t-1]
  }
}
plot(mixden.rescale(xsam.orig))
hist(xsam.orig, breaks=20)

xsam.rescale <- 0
for(t in 2:n.iter){
  prop <- rnorm(1,sd=1)+xsam.rescale[t-1]
  if(log(runif(1)) < mixden.rescale(prop)-mixden.rescale(xsam.rescale[t-1])){
    xsam.rescale[t] <- prop
  }else{
    xsam.rescale[t] <- xsam.orig[t-1]
  }
}
plot(mixden.rescale(xsam.rescale))
hist(xsam.rescale, breaks=200)
