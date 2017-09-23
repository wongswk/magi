library(parallel)
library(Rcpp)
sourceCpp("../src/wrapper.cpp")
source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

kerneltype <- "rbf" # this can be changed to "matern" for contrast

nobs <- 401
noise <- 0.05


VRtrue <- read.csv("../data/FN.csv")
pram.true <- list(
  abc=c(0.2,0.2,3),
  sigma=noise
)
fn.true <- VRtrue
fn.true$time <- seq(0,20,length=nrow(fn.true))
fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)

tvec.full <- fn.sim$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

phisigllikTest( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r, kerneltype)
fn <- function(par) -phisigllikTest( par, data.matrix(fn.sim[,1:2]), r, kerneltype)$value
gr <- function(par) -as.vector(phisigllikTest( par, data.matrix(fn.sim[,1:2]), r, kerneltype)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
marlikmap$par

curCovV <- calCov(marlikmap$par[1:2], kerneltype)
curCovR <- calCov(marlikmap$par[3:4], kerneltype)
cursigma <- marlikmap$par[5]

# marginal likelihood seems correct
phisigllik(marlikmap$par, data.matrix(fn.sim[,1:2]), TRUE, kerneltype)
phisigllikTest( marlikmap$par, data.matrix(fn.sim[,1:2]), r, kerneltype)
logliknoODE.mar(curCovV,
                curCovR, 
                cursigma,  
                data.matrix(fn.sim[,1:2]))
marlikmap$value
  
#' FIXME log full likelihood on derivative doesn't seem right, 
#' also prior smoothing part seems to have poor goodness of fit
xthetallik( data.matrix(fn.true[, 1:2]), 
            c(0.2, 0.2, 3),
            curCovV,
            curCovR, 
            cursigma,  
            data.matrix(fn.sim[,1:2]),
            F)
logliknoODE( data.matrix(fn.true[, 1:2]),
             curCovV,
             curCovR, 
             cursigma,  
             data.matrix(fn.sim[,1:2]))

loglik( data.matrix(fn.true[, 1:2]),
        c(0.2, 0.2, 3),
        curCovV,
        curCovR, 
        cursigma,  
        data.matrix(fn.sim[,1:2]))

#--- on sampled data, fitted derivative and true derivative are similar ----
#' the variance of derivative seems wrong

# V
draws <- MASS::mvrnorm(10, rep(0, nrow(fn.sim)), curCovV$C)
matplot(t(draws), type="l")
lines(data.matrix(fn.true[, 1]), lwd=2)

x <- draws[1,]
dxNum <- c(diff(head(x,2))/(fn.sim$time[2]-fn.sim$time[1]),
           (x[-(1:2)]-x[1:(length(x)-2)])/(fn.sim$time[3]-fn.sim$time[1]),
           diff(tail(x,2))/(fn.sim$time[2]-fn.sim$time[1]))
dxMean <- curCovV$mphi%*%x
plot(fn.sim$time, x, type="l", main="level")
plot(fn.sim$time, dxNum, type="l", main="derivative")
lines(fn.sim$time, dxMean, col=2)
plot(fn.sim$time, (dxNum-dxMean)/sqrt(diag(curCovV$Kphi)))

# R
draws <- MASS::mvrnorm(10, rep(0, nrow(fn.sim)), curCovR$C)
matplot(t(draws), type="l")
lines(data.matrix(fn.true[, 2]), lwd=2)

x <- draws[1,]
dxNum <- c(diff(head(x,2))/(fn.sim$time[2]-fn.sim$time[1]),
           (x[-(1:2)]-x[1:(length(x)-2)])/(fn.sim$time[3]-fn.sim$time[1]),
           diff(tail(x,2))/(fn.sim$time[2]-fn.sim$time[1]))
dxMean <- curCovR$mphi%*%x
plot(fn.sim$time, x, type="l", main="level")
plot(fn.sim$time, dxNum, type="l", main="derivative")
lines(fn.sim$time, dxMean, col=2)
plot(fn.sim$time, (dxNum-dxMean)/sqrt(diag(curCovR$Kphi)))

#--- on real data, fitted derivative and true derivative are quite different ----

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

# V
x <- fn.true$Vtrue
Ctrue <- curCovV$C - 1e-7 * diag( nrow(curCovV$C))
xfit <- Ctrue%*%solve(curCovV$C, fn.true$Vtrue)
dxNum <- fn.true$dVtrue
dxMean <- curCovV$mphi%*%x
plot(fn.true$time, x, type="l")
lines(fn.true$time, xfit, col=2)
plot(fn.sim$time, dxNum, type="l")
lines(fn.sim$time, dxMean, col=2)
plot((dxNum - dxMean)/sqrt(diag(curCovV$Kphi)), type="l")

# R
x <- fn.true$Rtrue
Ctrue <- curCovR$C - 1e-7 * diag( nrow(curCovR$C))
xfit <- Ctrue%*%solve(curCovR$C, fn.true$Rtrue)
dxNum <- fn.true$dRtrue
dxMean <- curCovR$mphi%*%x
plot(fn.true$time, x, type="l")
lines(fn.true$time, xfit, col=2)
plot(fn.sim$time, dxNum, type="l")
lines(fn.sim$time, dxMean, col=2)
plot((dxNum - dxMean)/sqrt(diag(curCovR$Kphi)), type="l")

