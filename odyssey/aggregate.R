source("../R/visualization.R")
source("../R/helper/utilities.r")
source("../R/helper/basic_hmc.R")
source("../R/HMC-functions.R")

load("../data/ody.RData")
args <- commandArgs(trailingOnly = TRUE)

temperature <- c(1,1,1)
lam <- 1/temperature # tuning parameter for weight on GP level fitting component

numparam <- nobs*2+3  # num HMC parameters
phi.ind <- seq(500,nrow(gpfit_ss$vphi), length.out=100)  ## which phi/sigma to use from initial fit
n.iter <- 300  # number of HMC iterations per phi
th.all <- matrix(NA,length(phi.ind),numparam)  # X and theta
phisig <- matrix(NA,length(phi.ind),5)   # phi and sigma

#th.all[1,] <-  c( startX, 1, 1, 1)
#th.all[1,] <-  c( startX, .2, .2, 3)
#th.all[1,] <- c( VRtrue[seq(1, 401, length = nobs),1], VRtrue[seq(1, 401, length = nobs),2], .2, .2, 3)
#phisig[1,] <- c( startphi, startsigma)
#phisig[1,] <- c( apply(gpfit_ss$vphi,2,median), apply(gpfit_ss$rphi,2,median), 0.5)

##### Reference values (truth)
ref.th <- c( VRtrue[seq(1, 401, length = nobs),1], VRtrue[seq(1, 401, length = nobs),2], .2, .2, 3)
bestCovV <- calCov( c(1.9840824, 1.1185157 ))
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
loglik( VRtrue[seq(1,401,length=nobs),], c(0.2,0.2,3), bestCovV, bestCovR, noise, fn.sim[,1:2], lambda=lam)
loglik( VRtrue[seq(1,401,length=nobs),], c(0.2,0.2,3), bestCovV, bestCovR, noise, fn.sim[,1:2], lambda=1)
loglik( VRtrue[seq(1,401,length=nobs),], c(0.2,0.2,3), bestCovV, bestCovR, noise*2, fn.sim[,1:2], lambda=4)

loglik( VRtrue[seq(1,401,length=nobs),], c(0.2,0.2,3), bestCovV, bestCovR, sigHigh, fn.sim[,1:2], lambda=lam)

## loglik at degenerate case (zero curve)
loglik(matrix(0,nrow=nobs,ncol=2),c(0,1,1),calCov(c(.1,10)), calCov(c(.1,10)), 0.25, fn.sim[,1:2], lambda=lam)
loglik(matrix(0,nrow=nobs,ncol=2),c(0,1,1),calCov(c(.1,10)), calCov(c(.1,10)), sigHigh, fn.sim[,1:2], lambda=lam)
loglik(matrix(0,nrow=nobs,ncol=2),c(0,1,1),calCov(c(.1,10)), calCov(c(.1,10)), noise, fn.sim[,1:2], lambda=lam)

## Bounds on phi and sigma
#lower_b <- c( 0, 0, 0, 0, sigLow )
#upper_b <- c( Inf, Inf, Inf, 10, sigHigh)
lower_b <- c( 0, 0, 0, 0, 0)
upper_b <- c( Inf, Inf, Inf, Inf, Inf)



#curllik <- xthU(th.all[1,], lambda = lam)

full_llik <- c()
lliklist <- c()
#lliklist[1] <- curllik
#full_llik[1] <- loglik( cbind(th.all[1,1:nobs],th.all[1,(nobs+1):(nobs*2)]), th.all[1,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)

#loglik( cbind(th.all[1,1:nobs],th.all[1,(nobs+1):(nobs*2)]), th.all[1,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
#accepts <- 0
#paccepts <- 0
#deltas <- c()


for (w in 1:length(phi.ind)){
  ret <- readRDS(paste0("/n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem", w, ".rds"))
  th.all[w,] <- ret[[1]]
  full_llik[w] <- ret[[2]]
  lliklist[w] <- ret[[3]]
  phisig[w,] <- c( gpfit_ss$vphi[phi.ind[w],], gpfit_ss$rphi[phi.ind[w],], gpfit_ss$sigma[phi.ind[w]])
}

id.best <- which.max(full_llik)
loglik( cbind(th.all[id.best,1:nobs],th.all[id.best,(nobs+1):(nobs*2)]), th.all[id.best,(nobs*2+1):(nobs*2+3)], calCov(phisig[id.best,1:2]),calCov(phisig[id.best,3:4]), phisig[id.best,5],  fn.sim[,1:2])

# pdf(file=paste0("R-HMC-output-",noise,".pdf"))
# par(mfrow=c(2,2))
# hist(th.all[501:5000,83], main="a")
# abline(v = 0.2, lwd=2, col="blue")
# hist(th.all[501:5000,84], main="b")
# abline(v = 0.2, lwd=2, col="blue")
# hist(th.all[501:5000,85], main="c")
# abline(v = 3, lwd=2, col="blue")
# hist(phisig[501:5000,5], main="sigma")
# abline(v = noise, lwd=2, col="blue")
# dev.off()

load("/Volumes/shihaoyang/Workspace/dynamicSystem/dynamic-systems/odyssey/ody_out.rda")

gpode <- list(abc=th.all[,(nobs*2+1):(nobs*2+3)],
              sigma=phisig[,5],
              rphi=phisig[,3:4],
              vphi=phisig[,1:2],
              rtrue=th.all[,(nobs+1):(nobs*2)],
              vtrue=th.all[,1:nobs],
              lp__=lliklist,
              lglik=full_llik)
gpode$fode <- sapply(1:length(phi.ind), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

plot.post.samples(paste0("../results/R-ode-",noise,".pdf"), fn.true, fn.sim, gpode, pram.true)

hist(gpode$abc[,1], breaks = 20)
hist(gpode$abc[,2], breaks = 20)
hist(gpode$abc[,3], breaks = 20)
