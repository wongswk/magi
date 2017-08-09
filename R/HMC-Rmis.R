### Required variables
### - gpfit_ss from the preliminary fit ignoring ODE
### - fn.sim with 41 rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)

### Load the variables above if necessary
#load("gpfit_ss.rda")
#load("fndata.rda")

source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

numparam <- 41*2+3  # num HMC parameters
n.iter <- 5000  # number of HMC iterations
th.all <- matrix(NA,n.iter,numparam)  # X and theta
phisig <- matrix(NA,n.iter,5)   # phi and sigma

#### Use gpfit_ss to set initial values
th.all[1,] <-  c( apply(gpfit_ss$vtrue, 2, mean), rep(0,41), 1, 1, 1)
phisig[1,] <- c( apply(gpfit_ss$vphi, 2, mean), 1, 1, mean(gpfit_ss$sigma))

##### Reference values (truth)
ref.th <- c( VRtrue[seq(1, 401, length = 41),1], VRtrue[seq(1, 401, length = 41),2], .2, .2, 3)
bestCovV <- calCov( c(1.9840824, 1.1185157 ))
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
loglikRmis( VRtrue[seq(1,401,length=41),], c(0.2,0.2,3), bestCovV, bestCovR, 0.1,  fn.sim[,1:2])

## loglik at degenerate case (zero curve)
loglikRmis(matrix(0,nrow=41,ncol=2),c(0,1,1),calCov(c(.1,10)), calCov(c(.1,10)), 0.25, fn.sim[,1:2])

## Bounds on phi and sigma
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, 0.25)


curCovV <- calCov(phisig[1,1:2])
curCovR <- calCov(phisig[1,3:4])
cursigma <- phisig[1,5]
curllik <- xthURmis(th.all[1,])

full_llik <- c()
lliklist <- c()
lliklist[1] <- curllik
full_llik[1] <- loglikRmis( cbind(th.all[1,1:41],th.all[1,42:82]), th.all[1,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2])
accepts <- 0
paccepts <- 0
#deltas <- c()


for (t in 2:n.iter) {
  
  if (t %% 100 == 0) { show(c(t, full_llik[t-1], accepts/t, paccepts/t)) }
  
  # Update X and theta
  foo <- basic_hmc(xthURmis, step=runif(1,0.001,0.002), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  lliklist[t] <- foo$lpr
  th.all[t,] <- foo$final
  #deltas[t] <- foo$delta
  accepts <- accepts + foo$acc
  
  # Update phi and sigma using random-walk M-H.
  oldCovV <- curCovV
  oldCovR <- curCovR
  old_ll <- loglikRmis( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2])
  ps_prop <- phisig[t-1,] + rnorm(5, 0, 0.03 * phisig[1,])
  if( min(ps_prop - lower_b) > 0 && min(upper_b - ps_prop) > 0) {  # check bounds
    propCovV <- calCov(ps_prop[1:2])
    propCovR <- calCov(ps_prop[3:4])
    prop_ll <- loglikRmis( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], propCovV, propCovR, ps_prop[5], fn.sim[,1:2])
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
  
  full_llik[t] <- loglikRmis( cbind(th.all[t,1:41],th.all[t,42:82]), th.all[t,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2])  
}  

## Best sampled
id.best <- which.max(full_llik)
loglikRmis( cbind(th.all[id.best,1:41],th.all[id.best,42:82]), th.all[id.best,83:85], calCov(phisig[id.best,1:2]),calCov(phisig[id.best,3:4]), phisig[id.best,5],  fn.sim[,1:2])

pdf(file="R-HMC-output-Rmis.pdf")
par(mfrow=c(2,2))
hist(th.all[501:5000,83], main="a")
abline(v = 0.2, lwd=2, col="blue")
hist(th.all[501:5000,84], main="b")
abline(v = 0.2, lwd=2, col="blue")
hist(th.all[501:5000,85], main="c")
abline(v = 3, lwd=2, col="blue")
hist(phisig[501:5000,5], main="sigma")
abline(v = 0.1, lwd=2, col="blue")

par(mfrow=c(1,2))
plot(tvec41, VRtrue[seq(1, 401, length = 41),1], type='l', main="V")
points(tvec41, th.all[5000,1:41], col="blue")

plot(tvec41, VRtrue[seq(1, 401, length = 41),2], type='l', main="R")
points(tvec41, th.all[5000,42:82], col="blue")
dev.off()


