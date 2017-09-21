library(parallel)
library(Rcpp)

nobs.candidates <- (2:14)^2+1
noise.candidates <- seq(0.05, 1.5, 0.05)

stratAtTruth <- TRUE

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
load((paste0("/n/kou_lab/shihaoyang/dynamic_sys/HMC-v4-fillat201/C-v4-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda")))

sourceCpp("../src/wrapper.cpp")
xth.formal[1,] <- c(fn.true$Vtrue, fn.true$Rtrue, pram.true$abc)

for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  #rstep <- rstep / 10
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      xth.formal[t-1,], rstep, 1000, T)
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
  
  if( t %% 100 == 0) show(c(t, mean(tail(accepts[1:t],100)), foo$final[(nall*2+1):(nall*2+3)]))
}



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
plot.post.samples(paste0("../results/postprocessed-HMC-v4-fixphi-noise",noise,"-nobs",nobs,"-startAtTruth.pdf"), fn.true, fn.sim, gpode, pram.true)
mean(accepts)
mean(stepLow.scaler)
saveRDS(gpode, paste0("../results/postprocessed-C-v4-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-startAtTruth.rds"))
