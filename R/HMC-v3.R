load("low_noise.RData")
Rcpp::sourceCpp('../src/wrapper.cpp')
source("visualization.R")
source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

phisigllikTest( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
fn <- function(par) -phisigllikTest( par, data.matrix(fn.sim[,1:2]), r)$value
gr <- function(par) -as.vector(phisigllikTest( par, data.matrix(fn.sim[,1:2]), r)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
marlikmap$par

curCovV <- calCov(marlikmap$par[1:2])
curCovR <- calCov(marlikmap$par[3:4])
cursigma <- marlikmap$par[5]

startVR <- rbind(getMeanCurve(fn.sim$time, fn.sim$Vtrue, fn.sim$time, 
                              t(marlikmap$par[1:2]), sigma.mat=matrix(cursigma)),
                 getMeanCurve(fn.sim$time, fn.sim$Rtrue, fn.sim$time, 
                              t(marlikmap$par[3:4]), sigma.mat=matrix(cursigma)))
startVR <- t(startVR)
nfold.pilot <- 5
nobs.pilot <- nobs%/%nfold.pilot
numparam.pilot <- nobs.pilot*2+3  # num HMC parameters
n.iter.pilot <- 1000
stepLow.scaler.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)
stepLow.scaler.pilot[1,] <-  0.01
accepts.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)
accepts.pilot[1,] <- 0
xth.pilot <- array(NA, dim=c(numparam.pilot, n.iter.pilot, nfold.pilot))
lliklist.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)

#### pilot HMC ####
# quickly navigate the initial value to right region
id.pilot <- matrix(1:(nfold.pilot*nobs.pilot), nobs.pilot, byrow = T)
it.pilot <- 1
for(it.pilot in 1:nfold.pilot){
  accepts <- c()
  stepLow <- rep(stepLow.scaler.pilot[1,it.pilot], nobs.pilot*2+3)
  
  fn.sim.pilot <- fn.sim[id.pilot[,it.pilot],]
  
  pilotSignedDist <- outer(fn.sim.pilot$time, t(fn.sim.pilot$time),'-')[,1,]
  
  xth.pilot[,1,it.pilot] <- c(startVR[id.pilot[,it.pilot],],rep(1,3))
  
  pilotCovV <- calCov2(marlikmap$par[1:2], abs(pilotSignedDist), -sign(pilotSignedDist))
  pilotCovR <- calCov2(marlikmap$par[3:4], abs(pilotSignedDist), -sign(pilotSignedDist))
  t <- 2
  for (t in 2:n.iter.pilot) {
    rstep <- runif(length(stepLow), stepLow, 2*stepLow)
    foo <- xthetaSample(data.matrix(fn.sim.pilot[,1:2]), pilotCovV, pilotCovR, cursigma, 
                        xth.pilot[,t-1,it.pilot], rstep, 20, T)
    xth.pilot[,t,it.pilot] <- foo$final
    accepts.pilot[t,it.pilot] <- foo$acc
    stepLow.scaler.pilot[t,it.pilot] <- tail(stepLow,1)
    if (mean(tail(accepts.pilot[1:t,it.pilot],100)) > 0.8) {
      stepLow <- stepLow * 1.01
    }else if (mean(tail(accepts.pilot[1:t,it.pilot],100)) < 0.5) {
      stepLow <- stepLow * .99
    }
    lliklist.pilot[t,it.pilot] <- foo$lpr
    if( t %% 100 == 0) show(c(t, mean(tail(accepts.pilot[1:t,it.pilot],100)), 
                              foo$final[(nobs.pilot*2+1):(nobs.pilot*2+3)]))
  }  
  cat("================\npilot iteration ", it.pilot, " finished\n================\n")
}

pdf("../results/pilot.pdf", height = 3*6, width = 5*6)
layout(matrix(1:15,3))
x <- apply(xth.pilot[-(1:(2*nobs.pilot)),,], c(1,3), plot)
dev.off()

#### setup parameters from pilot run ####
startTheta <- apply(xth.pilot[(nobs.pilot*2+1):(nobs.pilot*2+3),-(1:(dim(xth.pilot)[[2]]/2)),], 1, mean)
stepLow <- mean(stepLow.scaler.pilot[rbind(0,diff(accepts.pilot))==1])/nfold.pilot^2


#### formal running ####
numparam <- nobs*2+3
n.iter <- 5e3
xth.formal <- matrix(NA, numparam, n.iter)
xth.formal[,1] <- c(startVR,startTheta)
lliklist <- stepLow.scaler <- accepts <- rep(NA, n.iter)
accepts[1] <- 0
stepLow.scaler[1] <- stepLow
stepLow <- rep(stepLow, 2*nobs+3)
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[,t-1], rstep, 20, T)
  xth.formal[,t] <- foo$final
  accepts[t] <- foo$acc
  stepLow.scaler[t] <- tail(stepLow,1)
  
  if (t < n.iter/2) {
    if (mean(tail(accepts[1:t],100)) > 0.8) {
      stepLow <- stepLow * 1.01
    } else if (mean(tail(accepts[1:t],100)) < 0.5) {
      stepLow <- stepLow * .99
    }
  }
  lliklist[t] <- foo$lpr
  
  if( t %% 100 == 0) show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nobs*2+1):(nobs*2+3)]))
}
