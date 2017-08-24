setwd("~/Workspace/DynamicSys/dynamic-systems/test/")
library(Rcpp)
sourceCpp("../src/wrapper.cpp")
main()

out.all.c <- matrix(nrow=1e5, ncol=4)
out.all.c[1,] <- 1:4
for(i in 2:nrow(out.all.c)){
  out.normal <- hmc(out.all.c[i-1,], rep(0.05,4), 20, TRUE)
  out.all.c[i,] <- out.normal$final
}

hist(out.all.c[,4], breaks=100, probability = T)
plot.function(dnorm, from=-4,to=4,col=2,add=T)

hist(out.all.c[,3], breaks=100, probability = T)
plot.function(dnorm, from=-4,to=4,col=2,add=T)

hist(out.all.c[,2], breaks=100, probability = T)
plot.function(dnorm, from=-4,to=4,col=2,add=T)

hist(out.all.c[,1], breaks=100, probability = T)
plot.function(dnorm, from=-4,to=4,col=2,add=T)
