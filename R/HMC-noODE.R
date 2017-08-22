### Required variables
### - fn.sim with nobs rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)

### Load the variables above if necessary
#load("fndata.rda")
VRtrue <- read.csv("../data/FN.csv")

pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157)
)

nobs <- 50


library(parallel)
source("visualization.R")
source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

#noise level
noise <- 1.5
pram.true$sigma <- noise

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

#VRtrue[seq(1,401,length=nobs),] + rnorm(82, 0, noise) - VRtrue[seq(1,401,length=nobs),] 

numparam <- nobs*2  # num HMC parameters (no theta)
n.iter <- 5000  # number of HMC iterations
th.all <- matrix(NA,n.iter,numparam)  # X and theta
phisig <- matrix(NA,n.iter,5)   # phi and sigma

th.all[1,] <-  rep(0,nobs*2)
phisig[1,] <- rep(1,5)

##### Reference values (truth)
#ref.th <- c( VRtrue[seq(1, 401, length = nobs),1], VRtrue[seq(1, 401, length = nobs),2], .2, .2, 3)
bestCovV <- calCov( c(1.9840824, 1.1185157 ))
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
logliknoODE( VRtrue[seq(1,401,length=nobs),], bestCovV, bestCovR, 0.1,  fn.sim[,1:2])

## loglik at degenerate case (zero curve)
logliknoODE(matrix(0,nrow=nobs,ncol=2),calCov(c(.1,10)), calCov(c(.1,10)), 0.25, fn.sim[,1:2])

## Bounds on phi and sigma
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, Inf)


curCovV <- calCov(phisig[1,1:2])
curCovR <- calCov(phisig[1,3:4])
cursigma <- phisig[1,5]
curllik <- xthUnoODE(th.all[1,])

full_llik <- c()
lliklist <- c()
lliklist[1] <- curllik
full_llik[1] <- logliknoODE( cbind(th.all[1,1:nobs],th.all[1,(nobs+1):(nobs*2)]), curCovV, curCovR, cursigma,  fn.sim[,1:2])
accepts <- 0
paccepts <- 0
#deltas <- c()


gpmcmc <- mclapply(1:8, function(dummy.chain){
# th.all[1,] <-  rnorm(82)
# phisig[1,] <- exp(rnorm(5))
for (t in 2:n.iter) {
  
  if (t %% 100 == 0) { show(c(t, full_llik[t-1], accepts/t, paccepts/t)) }
  
  # Update X and theta
  foo <- basic_hmc(xthUnoODE, step=runif(1,0.002,0.004), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  lliklist[t] <- foo$lpr
  th.all[t,] <- foo$final
  #deltas[t] <- foo$delta
  accepts <- accepts + foo$acc
  
  # Update phi and sigma using random-walk M-H.
  oldCovV <- curCovV
  oldCovR <- curCovR
  old_ll <- logliknoODE( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2])
  ps_prop <- phisig[t-1,] + rnorm(5, 0, 0.05 * phisig[1,])
  if( min(ps_prop - lower_b) > 0 && min(upper_b - ps_prop) > 0) {  # check bounds
    propCovV <- calCov(ps_prop[1:2])
    propCovR <- calCov(ps_prop[3:4])
    prop_ll <- logliknoODE( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), propCovV, propCovR, ps_prop[5], fn.sim[,1:2])
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
  
  full_llik[t] <- logliknoODE( cbind(th.all[t,1:nobs],th.all[t,(nobs+1):(nobs*2)]), curCovV, curCovR, cursigma,  fn.sim[,1:2])  
}  
return(list(
  full_llik=full_llik,
  th.all=th.all,
  phisig=phisig,
  lliklist=lliklist
))
}, mc.cores = 8)

burnin <- 2500-1
full_llik <- do.call(c,lapply(gpmcmc, function(x) x$full_llik[-(1:burnin)]))
th.all <- do.call(rbind,lapply(gpmcmc, function(x) x$th.all[-(1:burnin),]))
phisig <- do.call(rbind,lapply(gpmcmc, function(x) x$phisig[-(1:burnin),]))
lliklist <- do.call(c,lapply(gpmcmc, function(x) x$lliklist[-(1:burnin)]))


## Best sampled
id.best <- which.max(full_llik)
logliknoODE( cbind(th.all[id.best,1:nobs],th.all[id.best,(nobs+1):(nobs*2)]), calCov(phisig[id.best,1:2]),calCov(phisig[id.best,3:4]), phisig[id.best,5],  fn.sim[,1:2])

# pdf(file="R-HMC-output.pdf")
# par(mfrow=c(2,2))
# hist(th.all[501:5000,83], main="a")
# abline(v = 0.2, lwd=2, col="blue")
# hist(th.all[501:5000,84], main="b")
# abline(v = 0.2, lwd=2, col="blue")
# hist(th.all[501:5000,85], main="c")
# abline(v = 3, lwd=2, col="blue")
# hist(phisig[501:5000,5], main="sigma")
# abline(v = 0.1, lwd=2, col="blue")
# dev.off()


show(quantile(phisig[2500:5000,5],c(0.001,0.999)))

startX <- apply(th.all[2500:5000,],2,mean)
startphi <- apply(phisig[2500:5000,1:4], 2, mean)
startsigma <- mean(phisig[2500:5000,5])
sigLow <- quantile(phisig[2500:5000,5], 0.001)
sigHigh <- quantile(phisig[2500:5000,5], 0.999)


gpfit <- list(sigma=phisig[,5],
              rphi=phisig[,3:4],
              vphi=phisig[,1:2],
              rtrue=th.all[,(nobs+1):(nobs*2)],
              vtrue=th.all[,1:nobs],
              lp__=lliklist,
              lglik=full_llik)



post.noODE <- summary.post.noODE(paste0("../results/R-GPfit-",noise,".pdf"), fn.true, fn.sim, gpfit, pram.true)

post.noODE$init.epost

#### check with STAN ####
gpfit.stan <- stan(file="../stan/gp-initialfit.stan",
              data=list(N=nrow(fn.sim),
                        robs=fn.sim$Rtrue,
                        vobs=fn.sim$Vtrue,
                        time=fn.sim$time),
              iter=600, chains=7, warmup = 200, cores=7)
gpfit_ss <- extract(gpfit.stan, permuted=TRUE)

post.noODE.stan <- summary.post.noODE(paste0("../results/STAN-noODE-noise-",noise,".pdf"), 
                                      fn.true, fn.sim, gpfit_ss, pram.true)

startX <- with(post.noODE.stan$init.epost, c(vtrue,rtrue))
startphi <- with(post.noODE.stan$init.epost, c(vphi,rphi))
startsigma <- post.noODE.stan$init.epost$sigma
sigLow <- exp(post.noODE.stan$gpfit.post["sigma","mean"]-4.5*post.noODE.stan$gpfit.post["sigma","sd"])
sigHigh <- exp(post.noODE.stan$gpfit.post["sigma","mean"]+4.5*post.noODE.stan$gpfit.post["sigma","sd"])

sigLow <- quantile(gpfit_ss$sigma, 0.001)

sigHigh <- quantile(gpfit_ss$sigma, 0.999)
