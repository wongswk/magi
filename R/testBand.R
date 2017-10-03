system("R CMD SHLIB ../src/bandmatlik.cpp -lRlapack -lRblas")
dyn.load("../src/bandmatlik.so")  # change to DLL if on Windows

source("../odyssey/HMC-v4-odyssey.R")

# Converts a matrix to banded form.
mat2band <- function(a, bandsize) {
  N <- nrow(a)
  A <- matrix(0,nrow = 2*bandsize+1, ncol = N)
  
  for (j in 1:N) {
    k <- bandsize + 1 - j
    for (i in max(1,j-bandsize):min(N,j+bandsize)) {
      A[k+i,j] <- a[i,j] 
    }
  }

  as.double(A)
}

# R wrapper C banded likelihood
xthetallik <- function(xtheta, Vmphi, VKinv, VCinv, Rmphi, RKinv, RCinv, bandsize, n, sigma,
                       yobs) {
  
  foo <- .C("xthetallik", xtheta = as.double(xtheta), 
     Vmphi = as.double(Vmphi), VKinv = as.double(VKinv), VCinv = as.double(VCinv),
     Rmphi = as.double(Rmphi), RKinv = as.double(RKinv), RCinv = as.double(RCinv),
     bandsize = as.integer(bandsize), nn = as.integer(n), sigma = as.double(sigma),
     yobs = as.double(yobs), ret = double(1), retgrad = double(length(xtheta)))
  ret <- foo$ret
  attr(ret,"grad") <- foo$retgrad
  
  ret
}

# Example call
bandsize <- 20
Vmphi <- mat2band(curCovV$mphi, bandsize)
VKinv <- mat2band(curCovV$Kinv, bandsize)
VCinv <- mat2band(curCovV$Cinv, bandsize)
Rmphi <- mat2band(curCovR$mphi, bandsize)
RKinv <- mat2band(curCovR$Kinv, bandsize)
RCinv <- mat2band(curCovR$Cinv, bandsize)
fn.sim.c <- fn.sim <- as.matrix(fn.sim)
fn.sim.c[is.nan(fn.sim.c)] <- -99999  # can't use NaN for missing values to pass to C.

approxCovV <- truncCovByEigen(curCovV, 
                              truncEigen(rev(curCovV$Ceigen1over), 0.85),
                              truncEigen(rev(curCovV$Keigen1over), 0.85))

approxCovR <- truncCovByEigen(curCovR, 
                              truncEigen(rev(curCovR$Ceigen1over), 0.85),
                              truncEigen(rev(curCovR$Keigen1over), 0.85))

out.appr <- xthetallik(xth.formal[1,], Vmphi, VKinv, VCinv, Rmphi, RKinv, RCinv, bandsize, 201, cursigma, fn.sim.c[,1:2])
out.orig <- xthetallikTest(fn.sim[,1:2], curCovV, curCovR, cursigma, xth.formal[1,])
out.lowr <- xthetallikTest(fn.sim[,1:2], approxCovV, approxCovR, cursigma, xth.formal[1,])  

delta <- 1e-8
gradNum <- c()
xthInit <- xth.formal[1,]
for(it in 1:length(xthInit)){
  xthInit1 <- xthInit
  xthInit1[it] <- xthInit1[it] + delta
  gradNum[it] <- 
    (xthetallik(xthInit1, Vmphi, VKinv, VCinv, Rmphi, RKinv, RCinv, bandsize, 201, cursigma, fn.sim.c[,1:2]) -
       xthetallik(xthInit, Vmphi, VKinv, VCinv, Rmphi, RKinv, RCinv, bandsize, 201, cursigma, fn.sim.c[,1:2]))/delta
}
x <- (gradNum - attr(out.appr,"grad"))/abs(attr(out.appr,"grad"))
plot(x)
summary(x)

(as.numeric(out.appr) - out.orig$value)/out.orig$value
(out.lowr$value - out.orig$value)/out.orig$value

microbenchmark::microbenchmark(
  xthetallik(xth.formal[1,], Vmphi, VKinv, VCinv, Rmphi, RKinv, RCinv, bandsize, 201, cursigma, fn.sim.c[,1:2]),
  xthetallikTest(fn.sim[,1:2], curCovV, curCovR, cursigma, xth.formal[1,]),
  xthetallikTest(fn.sim[,1:2], approxCovV, approxCovR, cursigma, xth.formal[1,])  
)
