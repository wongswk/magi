### Required variables
### - Run HMC-noODE.R first with chosen noise level to generate data and get preliminary GP fit
### - fn.sim with nobs rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)

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

w <- args

phisig[w,] <- c( gpfit_ss$vphi[phi.ind[w],], gpfit_ss$rphi[phi.ind[w],], gpfit_ss$sigma[phi.ind[w]])

curCovV <- calCov(phisig[w,1:2])
curCovR <- calCov(phisig[w,3:4])
cursigma <- phisig[w,5]

stepLow <- 0.00001
th.temp <- matrix(NA, n.iter, numparam)
th.temp[1,] <- c( gpfit_ss$vtrue[phi.ind[w],], gpfit_ss$rtrue[phi.ind[w],], 1, 1, 1)

accepts <- 0  

for (t in 2:n.iter) {
  
  # if (t %% 10 == 0) { show(c(t, full_llik[t-1], accepts/t, paccepts/t, stepLow)) }
  
  # Update X and theta
  #foo <- basic_hmc(xthU, step=runif(1,0.004,0.008), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  xthU.tempered <- function(q, grad) xthU(q, grad, lambda=lam)
  #foo <- basic_hmc(xthU.tempered, step=runif(1,0.001,0.002), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  foo <- basic_hmc(xthU.tempered, step=runif(1,stepLow,2*stepLow), nsteps= 20, initial=th.temp[t-1,], return.traj = T)
  #slliklist[t] <- foo$lpr
  th.temp[t,] <- foo$final
  #deltas[t] <- foo$delta
  accepts <- accepts + foo$acc
  if (t < n.iter/2) {
    if (accepts/t > 0.8) {
      stepLow <- stepLow * 1.01
    }
    if (accepts/t < 0.5) {
      stepLow <- stepLow * .99
    }
  }
}

th.all[w,] <- th.temp[t,]
full_llik[w] <- loglik( cbind(th.temp[t,1:nobs],th.temp[t,(nobs+1):(nobs*2)]), th.temp[t,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
lliklist[w] <- foo$lpr
# show(c(w, full_llik[w], accepts/t, th.all[w,(nobs*2+1):(nobs*2+3)]))
ret <- list(th.all[w,], full_llik[w], lliklist[w])
saveRDS(ret, paste0("/n/regal/kou_lab/shihaoyang/dynamic_sys/DynamicSystem", w, ".rds"))


