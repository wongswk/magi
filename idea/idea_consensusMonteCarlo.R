library(parallel)
library(Rcpp)

nobs <- 290
noise <- 0.05

rdaname <- paste0("../idea/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda")
load(rdaname)

sourceCpp("../src/wrapper.cpp")
source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

plot.post.samples(paste0("../idea/HMC-v3-fixphi-noise",noise,"-nobs",nobs,".pdf"), fn.true, fn.sim, gpode, pram.true)

system(paste0("open ../idea/HMC-v3-fixphi-noise",noise,"-nobs",nobs,".pdf"))

nfold.pilot <- ceiling(nobs/50)
id.pilot <- matrix(c(1:nobs, rep(NA, nfold.pilot-nobs%%nfold.pilot)), 
                   ncol=nfold.pilot, byrow = F)
nobs.pilot <- nrow(id.pilot)
numparam.pilot <- nobs.pilot*2+3  # num HMC parameters
n.iter.pilot <- 5000
burnin.pilot <- 1000
stepLow.scaler.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)
stepLow.scaler.pilot[1,] <-  0.001
apr.pilot <- accepts.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)
apr.pilot[1,] <- accepts.pilot[1,] <- 0
rstep.pilot <- array(NA, dim=c(numparam.pilot, n.iter.pilot, nfold.pilot))
rstep.pilot[,1,] <- stepLow.scaler.pilot[1,]
xth.pilot <- array(NA, dim=c(numparam.pilot, n.iter.pilot, nfold.pilot))
lliklist.pilot <- matrix(NA, n.iter.pilot, nfold.pilot)

#### pilot HMC ####
# quickly navigate the initial value to right region
it.pilot <- 1
for(it.pilot in 1:nfold.pilot){
  nobs.pilot <- length(na.omit(id.pilot[,it.pilot]))
  accepts <- c()
  stepLow <- rep(stepLow.scaler.pilot[1,it.pilot], nobs.pilot*2+3)
  
  fn.sim.pilot <- fn.sim[na.omit(id.pilot[,it.pilot]),]
  
  pilotSignedDist <- outer(fn.sim.pilot$time, t(fn.sim.pilot$time),'-')[,1,]
  
  # xth.pilot[,1,it.pilot] <- c(startVR[id.pilot[,it.pilot],],rep(1,3))
  xth.pilot[,1,it.pilot] <- c(startVR[id.pilot[,it.pilot],],c(0.2,0.2,3))
  
  pilotCovV <- calCov2(marlikmap$par[1:2], abs(pilotSignedDist), -sign(pilotSignedDist))
  pilotCovR <- calCov2(marlikmap$par[3:4], abs(pilotSignedDist), -sign(pilotSignedDist))
  t <- 2
  for (t in 2:n.iter.pilot) {
    rstep <- runif(length(stepLow), stepLow, 2*stepLow)
    foo <- xthetaSample(data.matrix(fn.sim.pilot[,1:2]), pilotCovV, pilotCovR, cursigma, 
                        na.omit(xth.pilot[,t-1,it.pilot]), rstep, 400, T)
    xth.pilot[,t,it.pilot][is.finite(xth.pilot[,t-1,it.pilot])] <- foo$final
    accepts.pilot[t,it.pilot] <- foo$acc
    apr.pilot[t,it.pilot] <- foo$apr
    stepLow.scaler.pilot[t,it.pilot] <- tail(stepLow,1)
    rstep.pilot[,t,it.pilot][is.finite(xth.pilot[,t-1,it.pilot])] <- rstep
    if(t <= burnin.pilot){
      if (mean(tail(accepts.pilot[1:t,it.pilot],100)) > 0.9) {
        stepLow <- stepLow * 1.01
      }else if (mean(tail(accepts.pilot[1:t,it.pilot],100)) < 0.6) {
        stepLow <- stepLow * .99
      }  
    }
    lliklist.pilot[t,it.pilot] <- foo$lpr
    if( t %% 100 == 0) show(c(t, mean(tail(accepts.pilot[1:t,it.pilot],100)), 
                              foo$final[(nobs.pilot*2+1):(nobs.pilot*2+3)]))
  }  
  cat("================\npilot iteration ", it.pilot, " finished\n================\n")
}
nobs.pilot <- nrow(id.pilot)



pdf(paste0("../idea/HMC-v3-pilot-noise",noise,"-nobs",nobs,".pdf"), height = 3*6, width = 5*nfold.pilot)
layout(matrix(1:(3*nfold.pilot),3))
x <- apply(xth.pilot[-(1:(2*nobs.pilot)),-(1:burnin.pilot),,drop=FALSE], c(1,3), plot)
dev.off()
system(paste0("open ../idea/HMC-v3-pilot-noise",noise,"-nobs",nobs,".pdf"))


theta.pilot <- xth.pilot[(dim(xth.pilot)[[1]]-2):dim(xth.pilot)[[1]],,,drop=FALSE]
var.pilot <- lapply(1:nfold.pilot, function(fd) var(t(theta.pilot[,,fd])))
precisionmat.pilot <- lapply(var.pilot, solve)
consensus <- lapply(1:nfold.pilot, function(fd) precisionmat.pilot[[fd]]%*%theta.pilot[,,fd])
consensus <- solve(Reduce("+", precisionmat.pilot), Reduce("+", consensus))

hist(consensus[1,-(1:burnin.pilot)])
hist(consensus[2,-(1:burnin.pilot)])
hist(consensus[3,-(1:burnin.pilot)])

pdf(paste0("../idea/HMC-v3-pilot-noise",noise,"-nobs",nobs,"-consensus.pdf"), height = 6, width = 6)
hist(consensus[1,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC a", probability = TRUE, 
     xlim = range(consensus[1,-(1:(ncol(consensus)/2))], gpode$abc[,1]))
hist(gpode$abc[,1], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[1], col=2)


hist(consensus[2,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC b", probability = TRUE,
     xlim = range(consensus[2,-(1:(ncol(consensus)/2))], gpode$abc[,2]))
hist(gpode$abc[,2], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[2], col=2)

hist(consensus[3,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC c", probability = TRUE,
     xlim = range(consensus[3,-(1:(ncol(consensus)/2))], gpode$abc[,3]))
hist(gpode$abc[,3], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[3], col=2)

dev.off()
