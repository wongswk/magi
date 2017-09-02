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

nfold.pilot <- 5
nobs.pilot <- nobs%/%nfold.pilot
numparam.pilot <- nobs.pilot*2+3  # num HMC parameters
n.iter.pilot <- 1000
stepLow <- c(rep(0.0002, nobs.pilot*2), rep(0.0002,3))
stepLow.pilot <- matrix(NA, nobs.pilot*2+3, nfold.pilot)
xth.pilot <- array(NA, dim=c(numparam.pilot, n.iter.pilot, nfold.pilot))
lliklist.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)

#### pilot HMC ####
# quickly navigate the initial value to right region
id.pilot <- matrix(1:(nfold.pilot*nobs.pilot), nobs.pilot, byrow = T)
it.pilot <- 1
for(it.pilot in 1:nfold.pilot){
  accepts <- c()
  stepLow <- c(rep(0.0002, nobs.pilot*2), rep(0.0002,3))
  xth.pilot[,1,it.pilot] <- c(startX[c(id.pilot[,it.pilot], id.pilot[,it.pilot]+nobs)], rep(1,3))
  fn.sim.pilot <- fn.sim[id.pilot[,it.pilot],]
  
  pilotSignedDist <- outer(fn.sim.pilot$time, t(fn.sim.pilot$time),'-')[,1,]
  
  fn.pilot <- function(par)
    -phisigllikTest( par, data.matrix(fn.sim.pilot[,1:2]), abs(pilotSignedDist))$value
  gr.pilot <- function(par)
    -as.vector(phisigllikTest( par, data.matrix(fn.sim.pilot[,1:2]), abs(pilotSignedDist))$grad)
  marlikmap.pilot <- optim(rep(1,5), fn.pilot, gr.pilot, method="L-BFGS-B", lower = 0.0001)
  # marlikmap.pilot <- marlikmap
  
  pilotCovV <- calCov2(marlikmap.pilot$par[1:2], abs(pilotSignedDist), sign(pilotSignedDist))
  pilotCovR <- calCov2(marlikmap.pilot$par[3:4], abs(pilotSignedDist), sign(pilotSignedDist))
  
  for (t in 2:n.iter.pilot) {
    rstep <- runif(length(stepLow), stepLow, 2*stepLow)
    foo <- xthetaSample(data.matrix(fn.sim.pilot[,1:2]), pilotCovV, pilotCovR, cursigma, 
                        xth.pilot[,t-1,it.pilot], rstep, 20, T)
    xth.pilot[,t,it.pilot] <- foo$final
    accepts <- c(accepts, foo$acc)
    if (mean(tail(accepts,100)) > 0.8) {
      stepLow <- stepLow * 1.01
    }else if (mean(tail(accepts,100)) < 0.5) {
      stepLow <- stepLow * .99
    }
    lliklist.pilot[t,it.pilot] <- foo$lpr
    if( t %% 100 == 0) show(c(t, mean(tail(accepts,100)), 
                              foo$final[(nobs.pilot*2+1):(nobs.pilot*2+3)]))
  }  
  stepLow.pilot[,it.pilot] <- stepLow
  cat("================\npilot iteration ", it.pilot, " finished\n================\n")
}

pdf("pilot.pdf", height = 3*6, width = 5*6)
layout(matrix(1:15,3))
x <- apply(xth.pilot[-(1:(2*nobs.pilot)),,], c(1,3), plot)
dev.off()

