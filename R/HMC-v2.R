### Required variables
### - Run HMC-noODE.R first with chosen noise level to generate data and get preliminary GP fit
### - fn.sim with nobs rows (noisy V and R in cols 1 & 2, using sigma = 0.1)
### - VRtrue with 401 rows (V and R true)

source("visualization.R")
source("helper/utilities.r")
source("helper/basic_hmc.R")
source("HMC-functions.R")

temperature <- c(1,1,1)
lam <- 1/temperature # tuning parameter for weight on GP level fitting component
if(!exists("gpfit_ss")) gpfit_ss <- gpfit

numparam <- nobs*2+3  # num HMC parameters
phi.ind <- seq(500,nrow(gpfit_ss$vphi), length.out=100)  ## which phi/sigma to use from initial fit
n.iter <- 1500  # number of HMC iterations per phi
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

for (w in 1:length(phi.ind)) {
  
  phisig[w,] <- c( gpfit_ss$vphi[phi.ind[w],], gpfit_ss$rphi[phi.ind[w],], gpfit_ss$sigma[phi.ind[w]])

  curCovV <- calCov(phisig[w,1:2])
  curCovR <- calCov(phisig[w,3:4])
  cursigma <- phisig[w,5]
  
  stepLow <- 0.001
  th.temp <- matrix(NA, n.iter, numparam)
  th.temp[1,] <- c( gpfit_ss$vtrue[phi.ind[w],], gpfit_ss$rtrue[phi.ind[w],], 1, 1, 1)

  accepts <- 0  
  
  for (t in 2:n.iter) {
    
    #if (t %% 100 == 0) { show(c(t, full_llik[t-1], accepts/t, paccepts/t, stepLow)) }
    
    # Update X and theta
    #foo <- basic_hmc(xthU, step=runif(1,0.004,0.008), nsteps= 20, initial=th.all[t-1,], return.traj = T)
    xthU.tempered <- function(q, grad) xthU(q, grad, lambda=lam)
    #foo <- basic_hmc(xthU.tempered, step=runif(1,0.001,0.002), nsteps= 20, initial=th.all[t-1,], return.traj = T)
    # foo <- basic_hmc(xthU.tempered, step=runif(1,stepLow,2*stepLow), nsteps= 20, initial=th.temp[t-1,], return.traj = T)
    foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                        th.temp[t-1,], rep(runif(1,stepLow,2*stepLow),ncol(th.temp)), 20, T)
    
    
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
    
    
    # # Update phi and sigma using random-walk M-H.
    # oldCovV <- curCovV
    # oldCovR <- curCovR
    # old_ll <- loglik( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), th.all[t,(nobs*2+1):(nobs*2+3)], oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2], lambda=lam) + logEB(oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2])
    # #old_ll <- loglik( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), th.all[t,(nobs*2+1):(nobs*2+3)], oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2], lambda=lam) 
    # 
    # ps_prop <- phisig[t-1,] + rnorm(5, 0, 0.05 * phisig[1,])
    # #ps_prop <- phisig[t-1,] + rnorm(5, 0, c(rep(0.00,4),0.25) * phisig[1,])
    # #phiind <-  sample(nrow(gpfit_ss$vphi),1)
    # #ps_prop[1:2] <- gpfit_ss$vphi[phiind,]
    # #ps_prop[3:4] <- gpfit_ss$rphi[phiind,]
    # if( min(ps_prop - lower_b) > 0 && min(upper_b - ps_prop) > 0) {  # check bounds
    #   propCovV <- calCov(ps_prop[1:2])
    #   propCovR <- calCov(ps_prop[3:4])
    #   prop_ll <- loglik( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), th.all[t,(nobs*2+1):(nobs*2+3)], propCovV, propCovR, ps_prop[5], fn.sim[,1:2], lambda=lam) + logEB(propCovV, propCovR, ps_prop[5], fn.sim[,1:2])
    #   #prop_ll <- loglik( cbind(th.all[t,1:nobs], th.all[t,(nobs+1):(nobs*2)]), th.all[t,(nobs*2+1):(nobs*2+3)], propCovV, propCovR, ps_prop[5], fn.sim[,1:2], lambda=lam)
    # } else {
    #   prop_ll <- -1e9  # reject if outside bounds
    # }
    # 
    # if (runif(1) < min(1,exp(prop_ll - old_ll))) {
    #   phisig[t,] <- ps_prop
    #   curCovV <- propCovV
    #   curCovR <- propCovR
    #   cursigma <- ps_prop[5]
    #   paccepts <- paccepts + 1
    # } else {
    #   phisig[t,] <- phisig[t-1,]
    # }
    # 
    #full_llik[t] <- loglik( cbind(th.all[t,1:nobs],th.all[t,(nobs+1):(nobs*2)]), th.all[t,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
  }
  
  th.all[w,] <- th.temp[t,]
  full_llik[w] <- loglik( cbind(th.temp[t,1:nobs],th.temp[t,(nobs+1):(nobs*2)]), th.temp[t,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, cursigma,  fn.sim[,1:2], lambda=lam)
  lliklist[w] <- foo$lpr
  show(c(w, full_llik[w], accepts/t, th.all[w,(nobs*2+1):(nobs*2+3)]))
}

## Best sampled
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


#### fixing phi sigma at marginal likelihood ####
phisigllikTest( c(1.9840824, 1.1185157, 0.9486433, 3.2682434, noise), data.matrix(fn.sim[,1:2]), r)
fn <- function(par) -phisigllikTest( par, data.matrix(fn.sim[,1:2]), r)$value
gr <- function(par) -as.vector(phisigllikTest( par, data.matrix(fn.sim[,1:2]), r)$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)
marlikmap$par

curCovV <- calCov(marlikmap$par[1:2])
curCovR <- calCov(marlikmap$par[3:4])
cursigma <- marlikmap$par[5]

n.iter <- 5000
stepLow <- c(rep(0.0001, nobs*2), rep(0.0001,3))
th.temp <- matrix(NA, n.iter, numparam)
th.temp[1,] <- c( colMeans(gpfit_ss$vtrue), colMeans(gpfit_ss$rtrue), 0.2, 0.2, 3)
#' initiating at 1,1,1 was not working.
#' reason is HMC sampler stuck in weird local mode
#' local model problem is more severe when observation is large and error is small
#' need something like parallel tempering to make the chain move away from local mode
#' anyway this is a sampler problem, not a method problem
full_llik <- c()
lliklist <- c()

accepts <- 0  


for (t in 2:n.iter) {
  foo <- xthetaSample(data.matrix(fn.sim[,1:2]), curCovV, curCovR, cursigma, 
                      th.temp[t-1,], stepLow, 20, T)
  th.temp[t,] <- foo$final
  accepts <- accepts + foo$acc
  if (t < n.iter/2) {
    if (accepts/t > 0.8) {
      stepLow <- stepLow * 1.01
    }
    if (accepts/t < 0.5) {
      stepLow <- stepLow * .99
    }
  }
  full_llik[t] <- loglik( cbind(th.temp[t,1:nobs],th.temp[t,(nobs+1):(nobs*2)]), 
                          th.temp[t,(nobs*2+1):(nobs*2+3)], curCovV, curCovR, 
                          cursigma,  fn.sim[,1:2], lambda=lam)
  lliklist[t] <- foo$lpr
  
  if( t %% 100 == 0) show(c(t, accepts/t, foo$final[(nobs*2+1):(nobs*2+3)]))
}



burnin <- 500

gpode <- list(abc=th.temp[-(1:burnin),(nobs*2+1):(nobs*2+3)],
              sigma=rep(marlikmap$par[5], n.iter-burnin),
              rphi=matrix(marlikmap$par[3:4], ncol=2,nrow=n.iter-burnin,byrow=T),
              vphi=matrix(marlikmap$par[1:2], ncol=2,nrow=n.iter-burnin,byrow=T),
              rtrue=th.temp[-(1:burnin),(nobs+1):(nobs*2)],
              vtrue=th.temp[-(1:burnin),1:nobs],
              lp__=lliklist[-(1:burnin)],
              lglik=full_llik[-(1:burnin)])
gpode$fode <- sapply(1:length(gpode$lp__), function(t) 
  with(gpode, fODE(abc[t,], cbind(vtrue[t,],rtrue[t,]))), simplify = "array")

fn.true$dVtrue = with(c(fn.true,pram.true), abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(c(fn.true,pram.true), -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

plot.post.samples(paste0("../results/R-ode-",noise,".pdf"), fn.true, fn.sim, gpode, pram.true)
