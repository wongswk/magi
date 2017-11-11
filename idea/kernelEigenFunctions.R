# convergence of eigen value and eigen vector ----------------------------------
maxT <- 5
tAll <- seq(0, maxT, length=41)
eigenList <- lapply(10:length(tAll), function(ndis){
  tDis <- seq(0, maxT, length=ndis)
  signedDist <- outer(tDis, tDis, '-')
  gpcov <- calCov(c(1,1), abs(signedDist), -sign(signedDist))
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

# fourier series for super-imposed sin curves ----------------------------------
rm(list=ls())
acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time

dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components

f   <- function(t,w) { 
  dc.component + 
    sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}

plot.function(Vectorize(function(t) f(t, 2*pi*f.0)),from = 0, to = 6)

w <- 2*pi*f.0
trajectory <- sapply(ts, function(t) f(t,w))
head(trajectory,n=30)
points(ts, trajectory)
X.k <- fft(trajectory)

P <- 1/f.0
x0 <- 0
ndis <- length(X.k)
fourierFft <- X.k / ndis

Kreconstruct2 <- function(ts) Re(sapply(ts, function(tj) 
  sum(fourierFft * exp(1i * (0:(ndis-1)) * 2*pi * (tj-x0) / P))))

plot.function(Kreconstruct2, from = 0, to = 6, add=TRUE, col=3, n=1e3)

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

plot.frequency.spectrum(fourierFft, xlimits=c(0,20))

get.trajectory <- function(X.k,ts,acq.freq) {
  
  N   <- length(ts)
  x.n <- rep(0,N)           # create vector to keep the trajectory
  ks  <- 0:(length(X.k)-1)
  
  for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
    x.n[n+1] <- sum(X.k * exp(1i*2*pi*ks*n/N)) / N
  }
  
  x.n * acq.freq 
}

x.n <- get.trajectory(X.k,ts,acq.freq) / acq.freq  # TODO: why the scaling?
# this is simply discrete inverse fourier transform that I already know
plot(ts,x.n, type="l"); abline(h=0,lty=3)
lines(ts,trajectory,col="red") # compare with original

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
  Kfunc(t) * exp(-1i*omega*t)
}

omegaCandidates <- seq(-6, 6, 0.01)

fourierTransKnumerical <- 
  sapply(omegaCandidates, function(omg)
    integrate( function(s) Re(fourierTransMatern(s, omg)), -10, 10 )$value)


plot(omegaCandidates, fourierTransKnumerical, type="l")

loss <- function(par) 
  sum((par/(2*maternDf/maternParm[2]^2+omegaCandidates^2)^(maternDf+0.5) - 
         fourierTransKnumerical)^2)

parms <- optim(c(1), loss, method = "L-BFGS-B", lower = c(0))
constParms <- parms$par 
constParms
fourierTransK <- function(omega) 
  constParms/(2*maternDf/maternParm[2]^2+omega^2)^(maternDf+0.5)

fourierTransKanalytical <- fourierTransK(omegaCandidates)

lines(omegaCandidates, fourierTransKanalytical, col=2)

# fourier series for matern kernel ------------------------------------------
maxT <- 5
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4)
P <- maxT*2

Napprox <- 9
Cn <- fourierTransK(2*pi*(0:Napprox)/P)/P

Kreconstruct <- function(ti) {
  Cn[1] + 2*sum(Cn[-1]*cos(2*pi*(1:Napprox)*ti/P))
}
Kreconstruct <- Vectorize(Kreconstruct)

plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4)
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)

plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2)
mycolor <- rev(grDevices::gray.colors(9, end = 0.7))
for(Napprox in 1:9){
  Cn <- fourierTransK(2*pi*(0:Napprox)/P)/P
  plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4, 
                col=mycolor[Napprox], add=TRUE)
}
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)

maxT <- maxT*3
plot.function(Kreconstruct, from = -maxT, to = maxT, n=1e4)
plot.function(Kfunc, from = -maxT, to = maxT, n=1e4, col=2, add=TRUE)
