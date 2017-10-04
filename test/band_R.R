system("R CMD SHLIB ../src/bandmatlik.cpp -lRlapack -lRblas")
dyn.load("../src/bandmatlik.so")  # change to DLL if on Windows
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

datainput <- scan("data_band.txt", sep = "\n", what = character())
datainput <- strsplit(datainput, "\t")
datainput <- lapply(datainput, function(x) as.numeric(na.omit(as.numeric(x))))

foo <- xthetallik(xtheta = datainput[[1]], 
          Vmphi = datainput[[2]], VKinv = datainput[[3]], VCinv = datainput[[4]],
          Rmphi = datainput[[5]], RKinv = datainput[[6]], RCinv = datainput[[7]],
          bandsize = datainput[[8]], n = datainput[[9]], sigma = datainput[[10]],
          yobs = datainput[[11]])

as.numeric(foo)+sum(attr(foo, "grad"))

system("R CMD SHLIB ../test/band.cpp")
dyn.load("../test/band.so")  # change to DLL if on Windo

.C("main")

system("R CMD SHLIB ../test/sharecode.cpp")
dyn.load("../test/sharecode.so")  # change to DLL if on Windo

.C("mainsharecode")
.C("simple")

