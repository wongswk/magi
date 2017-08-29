setwd("~/Workspace/DynamicSys/dynamic-systems/test/")
library(Rcpp)
sourceCpp("../src/paralleltempering.cpp")

#### gaussian example ####
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

#### bimodal gaussian example ####
ret <- main3()
dim(ret)

temperature <- c(1, 1.5, 2, 3, 4.5, 6, 8) 
layout(matrix(1:4, 2))
for(id in 1:length(temperature)){
  hist(ret[2,id,], probability = T)
  plot.function(function(x)
    (dnorm(x, -4, sd=temperature[id]) + dnorm(x, 4, sd=temperature[id]))/2, 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[3,id,], probability = T)
  plot.function(function(x)
    (dnorm(x, -4, sd=temperature[id]) + dnorm(x, 4, sd=temperature[id]))/2, 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[4,id,], probability = T)
  plot.function(function(x)
    (dnorm(x, -4, sd=temperature[id]) + dnorm(x, 4, sd=temperature[id]))/2, 
    from=-20,to=20,col=2,add=T, n=1001)
  hist(ret[5,id,], probability = T)
  plot.function(function(x)
    (dnorm(x, -4, sd=temperature[id]) + dnorm(x, 4, sd=temperature[id]))/2, 
    from=-20,to=20,col=2,add=T, n=1001)
}
