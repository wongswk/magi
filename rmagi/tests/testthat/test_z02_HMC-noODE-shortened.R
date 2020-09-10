testthat::context("test run HMC-noODE")
### Required variables
### - fn.sim with nobs rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)
library(magi)
VRtrue <- read.csv(system.file("testdata/FN.csv", package="magi"))

phitrue <- list(
  compact1 = c(2.618, 6.381, 0.152, 9.636),
  rbf = c(0.838, 0.307, 0.202, 0.653),
  matern = c(2.04, 1.313, 0.793, 3.101),
  periodicMatern = c(2.04, 1.313, 9, 0.793, 3.101, 9)
)

nobs <- 41

set.seed(Sys.time())
kerneltype <- sample(c("compact1","rbf","matern"),1)

pram.true <- list(abc=c(0.2, 0.2, 3),
                  phi=phitrue[[kerneltype]])

#noise level
noise <- 0.01
pram.true$sigma <- noise

fn.true <- VRtrue
fn.true$time <- seq(0,20,0.05)
fn.sim <- fn.true

set.seed(123)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]

tvec.nobs <- fn.sim$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

n.iter <- 50  # number of HMC iterations
phisig <- matrix(NA,n.iter,length(phitrue[[kerneltype]])+1)   # phi and sigma

fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r, kerneltype)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r, kerneltype)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
marlikmap$par
c(pram.true$phi, pram.true$sigma)
-marlikmap$value
loglikAtTruth <- phisigllikC( c(pram.true$phi, pram.true$sigma), data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r, kerneltype)$value
test_that("maximum likelihood should be higher than value at true parameter",{
  expect_gt(-marlikmap$value, loglikAtTruth)
})

phisig[1,] <- marlikmap$par

##### Reference values (truth)
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, Inf)

full_llik <- c()
full_llik[1] <- phisigllikC( phisig[1,], data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r, kerneltype)$value
accepts <- 0
paccepts <- 0
yobs <- data.matrix(fn.sim[,1:2])

if(kerneltype=="matern"){
  stepLow <- 0.01
}else if(kerneltype=="rbf"){
  stepLow <- 0.01
}else if(kerneltype=="compact1"){
  stepLow <- 0.01
}

for (t in 2:n.iter) {
  foo <- phisigSample(data.matrix(fn.sim[,1:2]), r, phisig[t-1,],
                      rep(runif(1,stepLow,2*stepLow),5), 200, T, kerneltype)
  phisig[t,] <- foo$final
  accepts <- accepts + foo$acc
  full_llik[t] <- foo$lpr
}  

burnin <- n.iter/2

## Best sampled
id.best <- which.max(full_llik)

startphi <- apply(phisig[-(1:burnin),1:4], 2, mean)
startsigma <- mean(phisig[-(1:burnin),5])
sigLow <- quantile(phisig[-(1:burnin),5], 0.001)
sigHigh <- quantile(phisig[-(1:burnin),5], 0.999)

gpfit <- list(sigma=phisig[,5],
              rphi=phisig[,3:4],
              vphi=phisig[,1:2],
              lp__=full_llik,
              lglik=full_llik)

plotx <- seq(0,20,0.1)
gpfit$vtrue <- getMeanCurve(fn.sim$time, fn.sim$Vtrue, plotx, 
                            gpfit$vphi, sigma.mat=gpfit$sigma, kerneltype)
gpfit$rtrue <- getMeanCurve(fn.sim$time, fn.sim$Rtrue, plotx, 
                            gpfit$rphi, sigma.mat=gpfit$sigma, kerneltype)

post.noODE <- magi:::summary.post.noODE(paste0("C-GPfit-",noise,"-",kerneltype,".pdf"),
                                 fn.true, fn.sim, gpfit, pram.true, plotx)

startX <- c(post.noODE$init.epost$vtrue, post.noODE$init.epost$rtrue)

logliknoODEOutEpost <- logliknoODE( cbind(post.noODE$init.epost$vtrue, post.noODE$init.epost$rtrue),
                                calCov(post.noODE$init.epost$vphi, r, signr, kerneltype=kerneltype),
                                calCov(post.noODE$init.epost$rphi, r, signr, kerneltype=kerneltype),
                                post.noODE$init.epost$sigma,
                                data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]))
logliknoODEOutTrue <- logliknoODE( data.matrix(fn.true[seq(1,nrow(fn.true), length=nobs), 1:2]),
                                  calCov(pram.true$phi[1:2], r, signr, kerneltype=kerneltype),
                                  calCov(pram.true$phi[3:4], r, signr, kerneltype=kerneltype),
                                  pram.true$sigma,
                                  data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]))

# FIXME the rbf seems to have very full small likelihood
# but variance of rbf kernel is definitely correct
logliknoODEOutEpost 
logliknoODEOutTrue

