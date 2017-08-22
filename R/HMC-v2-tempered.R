### Required variables
### - Run HMC-noODE.R first with chosen noise level to generate data and get preliminary GP fit
### - fn.sim with 41 rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)

source("visualization.R")
source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

temperature <- c(0.1025,1,1)
lam <- 1/temperature # tuning parameter for weight on GP level fitting component

numparam <- 41*2+3  # num HMC parameters
n.iter <- 5000  # number of HMC iterations
th.all <- matrix(NA,n.iter,numparam)  # X and theta
phisig <- matrix(NA,n.iter,5)   # phi and sigma

th.all[1,] <-  c( startX, 1, 1, 1)
phisig[1,] <- c( startphi, startsigma)

##### Reference values (truth)
ref.th <- c( VRtrue[seq(1, 401, length = 41),1], VRtrue[seq(1, 401, length = 41),2], .2, .2, 3)
bestCovV <- calCov( c(1.9840824, 1.1185157 ))
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
loglik( VRtrue[seq(1,401,length=41),], c(0.2,0.2,3), bestCovV, bestCovR, noise, fn.sim[,1:2], lambda=lam)


## loglik at degenerate case (zero curve)
loglik(matrix(0,nrow=41,ncol=2),c(0,1,1),calCov(c(.1,10)), calCov(c(.1,10)), 1.2, fn.sim[,1:2], lambda=lam)

## Bounds on phi and sigma
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, sigHigh)


curCovV <- calCov(phisig[1,1:2])
curCovR <- calCov(phisig[1,3:4])
cursigma <- phisig[1,5]
curllik <- xthU(th.all[1,], lambda = lam)

full_llik <- c()
lliklist <- c()
lliklist[1] <- curllik
full_llik[1] <- loglik( cbind(th.all[1,1:41],th.all[1,42:82]), th.all[1,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
accepts <- 0
paccepts <- 0
#deltas <- c()


for (t in 2:n.iter) {
  
  if (t %% 100 == 0) { show(c(t, full_llik[t-1], accepts/t, paccepts/t)) }
  
  # Update X and theta
  #foo <- basic_hmc(xthU, step=runif(1,0.004,0.008), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  xthU.tempered <- function(q, grad) xthU(q, grad, lambda=lam)
  foo <- basic_hmc(xthU.tempered, step=runif(1,0.001,0.002), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  lliklist[t] <- foo$lpr
  th.all[t,] <- foo$final
  #deltas[t] <- foo$delta
  accepts <- accepts + foo$acc
  
  # Update phi and sigma using random-walk M-H.
  oldCovV <- curCovV
  oldCovR <- curCovR
  old_ll <- loglik( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2], lambda=lam)
  ps_prop <- phisig[t-1,] + rnorm(5, 0, 0.05 * phisig[1,])
  if( min(ps_prop - lower_b) > 0 && min(upper_b - ps_prop) > 0) {  # check bounds
    propCovV <- calCov(ps_prop[1:2])
    propCovR <- calCov(ps_prop[3:4])
    prop_ll <- loglik( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], propCovV, propCovR, ps_prop[5], fn.sim[,1:2], lambda=lam)
  } else {
    prop_ll <- -1e9  # reject if outside bounds
  }
  
  if (runif(1) < min(1,exp(prop_ll - old_ll))) {
    phisig[t,] <- ps_prop
    curCovV <- propCovV
    curCovR <- propCovR
    cursigma <- ps_prop[5]
    paccepts <- paccepts + 1
  } else {
    phisig[t,] <- phisig[t-1,]
  }
  
  full_llik[t] <- loglik( cbind(th.all[t,1:41],th.all[t,42:82]), th.all[t,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
}  

## Best sampled
id.best <- which.max(full_llik)
loglik( cbind(th.all[id.best,1:41],th.all[id.best,42:82]), th.all[id.best,83:85], calCov(phisig[id.best,1:2]),calCov(phisig[id.best,3:4]), phisig[id.best,5],  fn.sim[,1:2])

pdf(file=paste0("R-HMC-output-",noise,".pdf"))
par(mfrow=c(2,2))
hist(th.all[501:5000,83], main="a")
abline(v = 0.2, lwd=2, col="blue")
hist(th.all[501:5000,84], main="b")
abline(v = 0.2, lwd=2, col="blue")
hist(th.all[501:5000,85], main="c")
abline(v = 3, lwd=2, col="blue")
hist(phisig[501:5000,5], main="sigma")
abline(v = noise, lwd=2, col="blue")
dev.off()


burnin <- 2500
gpode <- list(abc=th.all[-(1:burnin),83:85],
              sigma=phisig[-(1:burnin),5],
              rphi=phisig[-(1:burnin),3:4],
              vphi=phisig[-(1:burnin),1:2],
              rtrue=th.all[-(1:burnin),42:82],
              vtrue=th.all[-(1:burnin),1:41],
              lp__=lliklist[-(1:burnin)],
              lglik=full_llik[-(1:burnin)])
gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true <- VRtrue
fn.true$time <- seq(0,20,0.05)
fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

plot.post.samples(paste0("../results/R-ode-tempered-",noise,"-.pdf"), fn.true, fn.sim, gpode, pram.true)
