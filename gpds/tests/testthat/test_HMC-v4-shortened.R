testthat::context("test run HMC-v4")

library(testthat)
suppressMessages(library(gpds))

kerneltype <- "matern"

nobs <- 26
noise <- 0.05

VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=noise
)
fn.true <- VRtrue[seq(1,401,by=2),]   #### reference is 201 points
fn.true$time <- seq(0,20,0.1)
fn.sim <- fn.true

set.seed(123)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
tvec.full <- fn.sim$time
fn.sim.all <- fn.sim
fn.sim[-seq(1,nrow(fn.sim), length=nobs),] <- NaN
fn.sim.obs <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]
tvec.nobs <- fn.sim$time[seq(1,nrow(fn.sim), length=nobs)]

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)
fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
marlikmap$par

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=20, kerneltype=kerneltype)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=20, kerneltype=kerneltype)
cursigma <- marlikmap$par[5]

startVR <- rbind(getMeanCurve(fn.sim.obs$time, fn.sim.obs$Vtrue, fn.sim.obs$time, 
                              t(marlikmap$par[1:2]), sigma.mat=rep(cursigma,nrow(fn.sim.obs)), kerneltype),
                 getMeanCurve(fn.sim.obs$time, fn.sim.obs$Rtrue, fn.sim.obs$time, 
                              t(marlikmap$par[3:4]), sigma.mat=rep(cursigma,nrow(fn.sim.obs)), kerneltype))
startVR <- t(startVR)
nfold.pilot <- ceiling(nobs/40)
nobs.pilot <- nobs%/%nfold.pilot
numparam.pilot <- nobs.pilot*2+3  # num HMC parameters
n.iter.pilot <- 50
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
  
  fn.sim.pilot <- fn.sim.obs[id.pilot[,it.pilot],]
  
  pilotSignedDist <- outer(fn.sim.pilot$time, t(fn.sim.pilot$time),'-')[,1,]
  
  xth.pilot[,1,it.pilot] <- c(startVR[id.pilot[,it.pilot],],rep(1,3))
  
  pilotCovV <- calCov(marlikmap$par[1:2], abs(pilotSignedDist), -sign(pilotSignedDist), kerneltype=kerneltype)
  pilotCovR <- calCov(marlikmap$par[3:4], abs(pilotSignedDist), -sign(pilotSignedDist), kerneltype=kerneltype)
  
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
}

x <- apply(xth.pilot[-(1:(2*nobs.pilot)),,,drop=FALSE], c(1,3), plot)

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
nall <- nrow(fn.sim)
numparam <- nall*2+3
n.iter <- 50
xth.formal <- matrix(NA, n.iter, numparam)

# upsample with GP to get initial latent X's at non-observation points
obs.ind <- seq(1,nrow(fn.sim),length=nobs)
C.Vxx <- curCovV$C[obs.ind,obs.ind]
C.Vxmx <- curCovV$C[-obs.ind,obs.ind]
C.Vxmxm <- curCovV$C[-obs.ind,-obs.ind]
C.Rxx <- curCovR$C[obs.ind,obs.ind]
C.Rxmx <- curCovR$C[-obs.ind,obs.ind]
C.Rxmxm <- curCovR$C[-obs.ind,-obs.ind]
CV <- C.Vxmxm - C.Vxmx %*% solve(C.Vxx) %*% t(C.Vxmx)
CR <- C.Rxmxm - C.Rxmx %*% solve(C.Rxx) %*% t(C.Rxmx)
CV[upper.tri(CV)] <- t(CV)[upper.tri(CV)]
CR[upper.tri(CR)] <- t(CR)[upper.tri(CR)]

V_upsamp <- MASS::mvrnorm(1,C.Vxmx %*% solve(C.Vxx) %*% startVR[,1], CV)
R_upsamp <- MASS::mvrnorm(1,C.Rxmx %*% solve(C.Rxx) %*% startVR[,2], CR)

fullstartVR <- rep(NA, numparam-3)
fullstartVR[c(obs.ind,nall+obs.ind)] <- startVR
fullstartVR[-c(obs.ind,nall+obs.ind)] <- c(V_upsamp, R_upsamp)

xth.formal[1,] <- c(fullstartVR,startTheta)
xth.formal[1,] <- c(fn.true$Vtrue,fn.true$Rtrue, pram.true$abc)
lliklist <- stepLow.scaler <- accepts <- rep(NA, n.iter)
accepts[1] <- 0
stepLow.scaler[1] <- stepLow

fullscalefac <- rep(NA, numparam)
fullscalefac[c(obs.ind,nall+obs.ind,(2*nall+1):(2*nall+3))] <- scalefac   ## fill in nearest scales for HMC step size
for (ii in 2:numparam) { 
  if (is.na(fullscalefac[ii])) {
    fullscalefac[ii] <- fullscalefac[ii-1]
  }
}

stepLow <- stepLow.scaler[1]*fullscalefac * nobs/nall*0.4
burnin <- as.integer(n.iter*0.3)
n.iter.approx <- n.iter*0.9
for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  #rstep <- rstep / 10
  if(t < n.iter.approx){
    foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                        xth.formal[t-1,], rstep, 1000, T, loglikflag = "band")    
  }else{
    foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                        xth.formal[t-1,], rstep, 1000, T, loglikflag = "usual")
  }
  
  xth.formal[t,] <- foo$final
  accepts[t] <- foo$acc
  stepLow.scaler[t] <- mean(stepLow)
  
  if (t < burnin & t > 10) {
    if (mean(tail(accepts[1:t],100)) > 0.9) {
      stepLow <- stepLow * 1.005
    } else if (mean(tail(accepts[1:t],100)) < 0.6) {
      stepLow <- stepLow * .995
    }
  }
  lliklist[t] <- foo$lpr
  
   # show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nall*2+1):(nall*2+3)]))
}
# FIXME I see NaN if the chain was running long enough. Needs to debug


gpode <- list(abc=xth.formal[-(1:burnin), (nall*2+1):(nall*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=xth.formal[-(1:burnin), (nall+1):(nall*2)],
              vtrue=xth.formal[-(1:burnin), 1:nall],
              lp__=lliklist[-(1:burnin)],
              lglik=lliklist[-(1:burnin)])
gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

fn.sim$time <- fn.sim.all$time
gpds:::plot.post.samples(paste0("test-run-HMCv4-",kerneltype,".pdf"), 
                         fn.true, fn.sim, gpode, pram.true, 5)
