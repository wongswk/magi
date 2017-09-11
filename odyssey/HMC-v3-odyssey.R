library(parallel)
library(Rcpp)
sourceCpp("../src/wrapper.cpp")
source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

nobs.candidates <- (2:20)^2+1
noise.candidates <- seq(0.05, 1.5, 0.05)

nobs <- 21
noise <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if(length(args)>0){
  args <- as.numeric(args)
  noise <- noise.candidates[args%%length(noise.candidates)+1]
  args <- args%/%length(noise.candidates)
  nobs <- nobs.candidates[args%%length(nobs.candidates)+1]
  print(c(noise, nobs))
}

VRtrue <- read.csv("../data/FN.csv")
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=noise
)
fn.true <- VRtrue
fn.true$time <- seq(0,20,0.05)
fn.sim <- fn.true

fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]

tvec.nobs <- fn.sim$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

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
nfold.pilot <- ceiling(nobs/40)
nobs.pilot <- nobs%/%nfold.pilot
numparam.pilot <- nobs.pilot*2+3  # num HMC parameters
n.iter.pilot <- 1000
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
    apr.pilot[t,it.pilot] <- foo$apr
    stepLow.scaler.pilot[t,it.pilot] <- tail(stepLow,1)
    rstep.pilot[,t,it.pilot] <- rstep
    if (mean(tail(accepts.pilot[1:t,it.pilot],100)) > 0.9) {
      stepLow <- stepLow * 1.01
    }else if (mean(tail(accepts.pilot[1:t,it.pilot],100)) < 0.6) {
      stepLow <- stepLow * .99
    }
    lliklist.pilot[t,it.pilot] <- foo$lpr
    if( t %% 100 == 0) show(c(t, mean(tail(accepts.pilot[1:t,it.pilot],100)), 
                              foo$final[(nobs.pilot*2+1):(nobs.pilot*2+3)]))
  }  
  cat("================\npilot iteration ", it.pilot, " finished\n================\n")
}

pdf(paste0("../results/HMC-v3-pilot-noise",noise,"-nobs",nobs,".pdf"), height = 3*6, width = 5*nfold.pilot)
layout(matrix(1:(3*nfold.pilot),3))
x <- apply(xth.pilot[-(1:(2*nobs.pilot)),,,drop=FALSE], c(1,3), plot)
dev.off()

#### setup parameters from pilot run ####
startTheta <- apply(xth.pilot[(nobs.pilot*2+1):(nobs.pilot*2+3),-(1:(dim(xth.pilot)[[2]]/2)),], 1, mean)
stepLow <- mean(stepLow.scaler.pilot[rbind(0,diff(accepts.pilot))==1])/nfold.pilot^2
steptwopart <- apply(rstep.pilot, 1, function(x) sum(x*(apr.pilot>0.8))/sum(apr.pilot>0.8))
mean(steptwopart[1:(2*nobs.pilot)])
mean(tail(steptwopart,3))
scalefac <- apply(xth.pilot[,-(1:(dim(xth.pilot)[[2]]/2)),,drop=FALSE], c(1,3), sd)
scalefacV <- as.vector(t(scalefac[1:nobs.pilot,]))
scalefacV <- c(scalefacV, rep(mean(scalefacV), nobs-nobs.pilot*nfold.pilot))
scalefacR <- as.vector(t(scalefac[(nobs.pilot+1):(2*nobs.pilot),]))
scalefacR <- c(scalefacR, rep(mean(scalefacR), nobs-nobs.pilot*nfold.pilot))
scalefacTheta <- rowMeans(tail(scalefac,3))
scalefac <- c(scalefacV, scalefacR, scalefacTheta)
scalefac <- scalefac/mean(scalefac)

#### formal running ####
numparam <- nobs*2+3
n.iter <- 5e3
xth.formal <- matrix(NA, n.iter, numparam)
xth.formal[1,] <- c(startVR,startTheta)
lliklist <- stepLow.scaler <- accepts <- rep(NA, n.iter)
accepts[1] <- 0
stepLow.scaler[1] <- stepLow
stepLow <- stepLow.scaler[1]*scalefac
burnin <- as.integer(n.iter*0.3)
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, 20, T)
  xth.formal[t,] <- foo$final
  accepts[t] <- foo$acc
  stepLow.scaler[t] <- mean(stepLow)
  
  if (t < burnin) {
    if (mean(tail(accepts[1:t],100)) > 0.9) {
      stepLow <- stepLow * 1.01
    } else if (mean(tail(accepts[1:t],100)) < 0.6) {
      stepLow <- stepLow * .99
    }
  }
  lliklist[t] <- foo$lpr
  
  if( t %% 100 == 0) show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nobs*2+1):(nobs*2+3)]))
}



gpode <- list(abc=xth.formal[-(1:burnin), (nobs*2+1):(nobs*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=xth.formal[-(1:burnin), (nobs+1):(nobs*2)],
              vtrue=xth.formal[-(1:burnin), 1:nobs],
              lp__=lliklist[-(1:burnin)],
              lglik=lliklist[-(1:burnin)])
gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

plot.post.samples(paste0("../results/HMC-v3-fixphi-noise",noise,"-nobs",nobs,".pdf"), fn.true, fn.sim, gpode, pram.true)
mean(accepts)
mean(stepLow.scaler)
save.image(paste0("../results/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda"))