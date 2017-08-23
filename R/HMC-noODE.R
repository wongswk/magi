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

nobs <- 201


library(parallel)
source("visualization.R")
source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

#noise level
noise <- 0.5
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

n.iter <- 100  # number of HMC iterations
phisig <- matrix(NA,n.iter,5)   # phi and sigma

phisig[1,] <- rep(1,5)

##### Reference values (truth)
bestCovV <- calCov( c( 1.9840824, 1.1185157) )
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
logliknoODE.mar( bestCovV, bestCovR, noise, fn.sim[,1:2])
phisigllik( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), fn.sim[,1:2])

## loglik at degenerate case (zero curve)
logliknoODE.mar(calCov(c(.1,10)), calCov(c(.1,10)), 1.25, fn.sim[,1:2])
phisigllik( c(.1,10,.1,10, 1.25), fn.sim[,1:2])

## Bounds on phi and sigma
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, Inf)

full_llik <- c()
full_llik[1] <- phisigllik( phisig[1,], fn.sim[,1:2])
accepts <- 0
paccepts <- 0
yobs <- data.matrix(fn.sim[,1:2])
phisigU <- function(phisigval, grad = F) phisigllik(phisigval, y = yobs, grad = grad)  

stepLow <- 0.01

gpmcmc <- mclapply(1:8, function(dummy.chain){
  
  for (t in 2:n.iter) {
    
    # if (t %% 10 == 0) { cat(c(t, full_llik[t-1], accepts/t), "\n") }

    foo <- basic_hmc(phisigU, step=runif(1,stepLow,2*stepLow), nsteps= 20, initial=phisig[t-1,], return.traj = T)
    phisig[t,] <- foo$final
    accepts <- accepts + foo$acc
    if (t < n.iter/2) {
      if (accepts/t > 0.8) stepLow <- stepLow * 1.01
      if (accepts/t < 0.5) stepLow <- stepLow * .99
    }
    full_llik[t] <- foo$lpr
  }  
  
  return(list(
    full_llik=full_llik,
    phisig=phisig
  ))
}, mc.cores = 8)

burnin <- 20
full_llik <- do.call(c,lapply(gpmcmc, function(x) x$full_llik[-(1:burnin)]))
phisig <- do.call(rbind,lapply(gpmcmc, function(x) x$phisig[-(1:burnin),]))


## Best sampled
id.best <- which.max(full_llik)
phisigllik( phisig[id.best,],  fn.sim[,1:2])

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


show(quantile(phisig[-(1:burnin),5],c(0.001,0.999)))

startphi <- apply(phisig[-(1:burnin),1:4], 2, mean)
startsigma <- mean(phisig[-(1:burnin),5])
sigLow <- quantile(phisig[-(1:burnin),5], 0.001)
sigHigh <- quantile(phisig[-(1:burnin),5], 0.999)


gpfit <- list(sigma=phisig[,5],
              rphi=phisig[,3:4],
              vphi=phisig[,1:2],
              # rtrue=th.all[,(nobs+1):(nobs*2)],
              # vtrue=th.all[,1:nobs],
              lp__=full_llik,
              lglik=full_llik)


gpfit$vtrue <- getMeanCurve(fn.sim$time, fn.sim$Vtrue, fn.sim$time, 
                            gpfit$vphi, sigma.mat=gpfit$sigma)
gpfit$rtrue <- getMeanCurve(fn.sim$time, fn.sim$Rtrue, fn.sim$time, 
                            gpfit$rphi, sigma.mat=gpfit$sigma)

startX <- colMeans(cbind(gpfit$vtrue, gpfit$rtrue))

post.noODE <- summary.post.noODE(paste0("../results/R-GPfit-",noise,".pdf"), fn.true, fn.sim, gpfit, pram.true)

post.noODE$init.epost

#### check with STAN ####
gpfit.stan <- stan(file="../stan/gp-initialfit-mar.stan",
                   data=list(N=nrow(fn.sim),
                             robs=fn.sim$Rtrue,
                             vobs=fn.sim$Vtrue,
                             time=fn.sim$time),
                   iter=600, chains=7, warmup = 200, cores=7)
gpfit_ss <- extract(gpfit.stan, permuted=TRUE)

gpfit_ss$vtrue <- getMeanCurve(fn.sim$time, fn.sim$Vtrue, fn.sim$time, 
                               gpfit_ss$vphi, sigma.mat=gpfit_ss$sigma)
gpfit_ss$rtrue <- getMeanCurve(fn.sim$time, fn.sim$Rtrue, fn.sim$time, 
                               gpfit_ss$rphi, sigma.mat=gpfit_ss$sigma)

post.noODE.stan <- summary.post.noODE(paste0("../results/STAN-noODE-noise-",noise,".pdf"), 
                                      fn.true, fn.sim, gpfit_ss, pram.true)

startX <- with(post.noODE.stan$init.epost, c(vtrue,rtrue))
startphi <- with(post.noODE.stan$init.epost, c(vphi,rphi))
startsigma <- post.noODE.stan$init.epost$sigma
sigLow <- exp(post.noODE.stan$gpfit.post["sigma","mean"]-4.5*post.noODE.stan$gpfit.post["sigma","sd"])
sigHigh <- exp(post.noODE.stan$gpfit.post["sigma","mean"]+4.5*post.noODE.stan$gpfit.post["sigma","sd"])

save.image(paste0("../results/STAN-noODE-noise-",noise,".RData"))

sigLow <- quantile(gpfit_ss$sigma, 0.001)

sigHigh <- quantile(gpfit_ss$sigma, 0.999)
