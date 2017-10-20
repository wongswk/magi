suppressMessages(library(gpds))

set.seed(Sys.time())
kerneltype <- "generalMatern"

nobs <- 41
noise <- 0.1

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

myseed <- sample(1e10, 1) # or use nano second system time
myseed <- 123
set.seed(myseed)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
tvec.full <- fn.sim$time
fn.sim.all <- fn.sim
fn.sim[-seq(1,nrow(fn.sim), length=nobs),] <- NaN
fn.sim.obs <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]
tvec.nobs <- fn.sim$time[seq(1,nrow(fn.sim), length=nobs)]

config <- list(
  nobs = nobs,
  noise = noise,
  kernel = kerneltype,
  seed = myseed,
  npostplot = 5
)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

phisigllikC( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)
fn <- function(par) -phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r.nobs, kerneltype)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
cursigma <- marlikmap$par[5]


n.iter <- 1000 # number of HMC iterations
phisig <- matrix(NA,n.iter,5)   # phi and sigma

phisig[1,] <- marlikmap$par

##### Reference values (truth)
lower_b <- c( 0, 0, 0, 0, 0 )
upper_b <- c( Inf, Inf, Inf, Inf, Inf)

full_llik <- c()
full_llik[1] <- phisigllikC( phisig[1,], data.matrix(fn.sim[!is.nan(fn.sim[,1]),1:2]), r, kerneltype)$value
accepts <- 0
paccepts <- 0
yobs <- data.matrix(fn.sim[,1:2])

stepLow <- c( rep(0.1,4),0.01)
accepts <- 0

for (t in 2:n.iter) {
  rstep <- runif(length(stepLow), stepLow, 2*stepLow)
  foo <- phisigSample(data.matrix(fn.sim.obs[,1:2]), r, phisig[t-1,],
                      rstep, 10, T, kerneltype)
  phisig[t,] <- foo$final
  accepts <- accepts + foo$acc
  full_llik[t] <- foo$lpr
}  

save(fn.sim,fn.sim.all,fn.sim.obs,phisig, VRtrue,tvec.full, tvec.nobs, file="varyphisig.rda")

