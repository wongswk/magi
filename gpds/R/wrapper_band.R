#' Converts a matrix to banded form
#' 
#' @export
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

#' R wrapper C banded likelihood
#' @export
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

