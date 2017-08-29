setwd("~/Workspace/DynamicSys/dynamic-systems/test/")
library(Rcpp)
sourceCpp("../src/paralleltempering.cpp")
ret <- main2()
dim(ret)

temperature <- 8:1
layout(matrix(1:4, 2))
for(id in 1:8){
  hist(ret[2,id,], probability = T)
  plot.function(function(x)
    dnorm(x, sd=sqrt(temperature[id])), 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[3,id,], probability = T)
  plot.function(function(x)
    dnorm(x, sd=sqrt(temperature[id])), 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[4,id,], probability = T)
  plot.function(function(x)
    dnorm(x, sd=sqrt(temperature[id])), 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[5,id,], probability = T)
  plot.function(function(x)
    dnorm(x, sd=sqrt(temperature[id])), 
    from=-20,to=20,col=2,add=T, n=1001)
}

