# convergence of eigen value and eigen vector ----------------------------------
maxT <- 5
tAll <- seq(0, maxT, length=41)
eigenList <- lapply(10:length(tAll), function(ndis){
  tDis <- seq(0, maxT, length=ndis)
  signedDist <- outer(tDis, tDis, '-')
  gpcov <- calCov(c(1,1), abs(signedDist), -sign(signedDist), kerneltype = "matern")
  outX <- apply(gpcov$CeigenVec[,1:10], 2, function(y) {
    y <- y * sqrt(ndis / maxT)
    outy <- approx(tDis, y, tAll)$y
    outy <- outy*sign(outy[1])
    outy
  })
  list(head(1/gpcov$Ceigen1over / ndis, 10), outX)
})
matplot(t(sapply(eigenList, function(x) x[[1]])), type="l")
eigenFun <- sapply(eigenList, function(x) x[[2]], simplify = "array")

mycolor <- rev(gray.colors(dim(eigenFun)[3]))
matplot(eigenFun[,1,], type="l", lty=1, col=mycolor)
matplot(eigenFun[,2,], type="l", lty=1, col=mycolor)
matplot(eigenFun[,3,], type="l", lty=1, col=mycolor)

# fourier series for sawtooth wave ------------------------------------------
Kfunc <- function(r) r/pi # defined on [-pi, pi]
P <- 2*pi
x0 <- -pi
maxT <- pi
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4)

ndis <- 40
tdis <- seq(-maxT, maxT, length=ndis+1)[-(ndis+1)]

xvals <- Kfunc(tdis)
fourierFft <- fft(xvals)/ndis

testthat::expect_equal(Re(sapply(0:(ndis-1), function(n)
  sum(fourierFft * exp(1i * 2*pi * 0:(ndis-1) * n / ndis)))), xvals)

testthat::expect_equal(Re(sapply(tdis, function(tj) 
  sum(fourierFft * exp(1i * (0:(ndis-1)) * 2*pi * (tj-x0) / P)))), xvals)

fourierCoef <- list(c0 = fourierFft[1])
fourierCoef$cPos <- 0.5 * fourierFft[-1] * exp(-1i*2*pi*(1:(ndis-1))*x0/P)
fourierCoef$cNeg <- Conj(fourierCoef$cPos)

Kreconstruct <- function(ti) {
  Re(fourierCoef$c0 + sum(fourierCoef$cPos * exp(1i * 2*pi * (1:length(fourierCoef$cPos)) * ti / P)) +
    sum(fourierCoef$cNeg * exp(-1i * 2*pi * (1:length(fourierCoef$cNeg)) * ti / P)))
}
Kreconstruct <- Vectorize(Kreconstruct)

Kreconstruct2 <- function(ts) Re(sapply(ts, function(tj) 
  sum(fourierFft * exp(1i * (0:(ndis-1)) * 2*pi * (tj-x0) / P))))

testthat::expect_equal(Kreconstruct(tdis), xvals)
testthat::expect_equal(Kreconstruct2(tdis), xvals)

fourierCoef <- list(
  c0 = 0,
  cPos = -0.5*1i*2*(-1)^(2:ndis)/(pi*(1:(ndis-1))),
  cNeg = 0.5*1i*2*(-1)^(2:ndis)/(pi*(1:(ndis-1)))
)

plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)
plot.function(Kreconstruct2, from = -maxT, to = maxT, n=1e4, col=3, add=TRUE)
# use discrete fft is not going to get back true fourier coefficients

# fourier transform for matern kernel ------------------------------------------
rm(list=ls())
maternDf <- 2.5
maternParm <- c(1,1)
Kfunc <- function(r) calCovGeneralMatern(maternParm, abs(r), NULL, complexity = 0, df=maternDf)$C


fourierTransMatern <- function(t, omega) {
  Kfunc(t) * exp(-1i*2*pi*omega*t)
}

omegaCandidates <- seq(-6, 6, 0.01)

fourierTransKnumerical <- 
  sapply(omegaCandidates, function(omg)
    integrate( function(s) Re(fourierTransMatern(s, omg)), -10, 10 )$value)

plot(omegaCandidates, fourierTransKnumerical, type="l")

constParms <- maternParm[1]*(2*maternDf/maternParm[2]^2)^(maternDf+0.5)*2*pi*maternParm[2]/
  (beta(0.5, maternDf)*sqrt(2*maternDf))
fourierTransK <- function(omega) 
  constParms/(2*maternDf/maternParm[2]^2+4*pi^2*omega^2)^(maternDf+0.5)

fourierTransKanalytical <- fourierTransK(omegaCandidates)

lines(omegaCandidates, fourierTransKanalytical, col=2)

testthat::expect_equal(integrate( Kfunc, -10, 10 )$value, fourierTransK(0), tol=1e-6)
testthat::expect_equal(integrate( fourierTransK, -10, 10 )$value, Kfunc(0), tol=1e-6)

# fourier series for matern kernel ------------------------------------------
maxT <- 5
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4)
P <- maxT*2

Napprox <- 9
Cn <- fourierTransK((0:Napprox)/P)/P

Kreconstruct <- function(ti) {
  Cn[1] + 2*sum(Cn[-1]*cos(2*pi*(1:Napprox)*ti/P))
}
Kreconstruct <- Vectorize(Kreconstruct)

plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4)
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)

plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2)
mycolor <- rev(grDevices::gray.colors(9, end = 0.7))
for(Napprox in 1:9){
  Cn <- fourierTransK((0:Napprox)/P)/P
  plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4, 
                col=mycolor[Napprox], add=TRUE)
}
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)

maxTexpand <- maxT*3
plot.function(Kreconstruct, from = -maxTexpand, to = maxTexpand, n=1e4)
plot.function(Kfunc, from = -maxTexpand, to = maxTexpand, n=1e4, col=2, add=TRUE)

plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4, ylim=c(-0.5,1))
for(n in 0:Napprox){
  plot.function(function(ti) (1+(n>0))*Cn[1+n]*cos(2*pi*n*ti/P),
    from = -maxT, to = maxT, n=1e4, col=3, add=T)
}

# eigen value/function from fourier series for matern kernel ------------------
# readings from "Hilbert Space Methods for Reduced-Rank Gaussian Process Regression"
L <- 1
phi <- function(x, j){
  sin(pi*j*(x+L)/(2*L))/sqrt(L)
}
plot.function(function(x) phi(x, 1), from = -L, to = L, n=1e3)  
plot.function(function(x) phi(x, 3), from = -L, to = L, n=1e3)  
