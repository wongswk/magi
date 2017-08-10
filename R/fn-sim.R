#### start of code ####
source("R/func.R")
source("R/visualization.R")
library(rstan)
fn.true <- read.csv("data/FN.csv")
fn.true$time <- seq(0,20,0.05)
matplot(fn.true$time, data.matrix(fn.true[,-3]), type="l", lty=1)

abc = c(0.2, 0.2, 3)

read.fn.sim <- TRUE
# set.seed(123)

fn.true$dVtrue = with(fn.true, abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(fn.true, -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=0.1)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=41),]

if(read.fn.sim){
  fn.sim <- read.csv("data/fn.sim.csv")[,-1]
}

matplot(fn.sim$time, data.matrix(fn.sim[,-3]), type="l", lty=1)

init <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=0.1
)

tvec41 <- fn.sim$time
foo <- outer(tvec41, t(tvec41),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)


init$reta <- solve(t(chol(calCov(init$rphi, r)$C)), fn.true[seq(1,401,length=41),c("Rtrue")])
init$veta <- solve(t(chol(calCov(init$vphi, r)$C)), fn.true[seq(1,401,length=41),c("Vtrue")])

stopifnot(abs(getX(r, t(init$rphi), t(init$reta)) - fn.true[seq(1,401,length=41),c("Rtrue")]) < 1e-14)

#### investigate true value ####
# not working on need to fix
hyperparm0 <- c(rep(0,5),rep(5,5))
hyperveta0 <- hyperreta0 <- cbind(rep(0,nrow(fn.sim)), rep(1,nrow(fn.sim)))

gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time,
                           hyperparm=hyperparm0,
                           hyperreta=hyperreta0, 
                           hyperveta = hyperveta0,
                           ubsigma = 1000),
                 iter=2, chains=1, init=list(init), warmup = 0)


gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)
gpsmooth_ss$lp__

# init has lp value 68.00089
# simulation around -180

stopifnot(abs(gpsmooth_ss$vtrue[1,] - fn.true[seq(1,401,length=41),c("Vtrue")])<1e-10)
stopifnot(abs(gpsmooth_ss$rtrue[1,] - fn.true[seq(1,401,length=41),c("Rtrue")])<1e-10)

covR.init <- calCov(init$rphi, r)
covV.init <- calCov(init$vphi, r)

Rpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$rphi, eta.mat = gpsmooth_ss$reta)
Vpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$vphi, eta.mat = gpsmooth_ss$veta)
dVdRpostsample <- getdVdR(abc.mat = gpsmooth_ss$abc, rtrue.mat = gpsmooth_ss$rtrue, vtrue.mat = gpsmooth_ss$vtrue)

#' checking value output at init is the same as true value
stopifnot(abs(gpsmooth_ss$C_rphi[1,,] - covR.init$C)<1e-12)
stopifnot(abs(gpsmooth_ss$L_C_rphi[1,,] - t(chol(covR.init$C))) < 1e-12)
stopifnot(abs(gpsmooth_ss$K_rphi[1,,] - covR.init$Kphi) < 1e-10)
stopifnot(abs(gpsmooth_ss$dC_rphi[1,,] - covR.init$Cprime) < 1e-12)
stopifnot(abs(gpsmooth_ss$ddC_rphi[1,,] - covR.init$Cdoubleprime) < 1e-12)
stopifnot(abs(gpsmooth_ss$m_rphi_rtrue[1,] - covR.init$mphi%*%fn.true[seq(1,401,length=41),c("Rtrue")]) < 1e-11)
stopifnot(abs(gpsmooth_ss$drobs[1,] - fn.true[seq(1,401,length=41),c("dRtrue")])<1e-10)
stopifnot(abs(Rpostsample - gpsmooth_ss$rtrue) < 1e-10)
stopifnot(abs(Vpostsample - gpsmooth_ss$vtrue) < 1e-10)
stopifnot(abs(dVdRpostsample[,,"drobs"] - gpsmooth_ss$drobs) < 1e-10)
stopifnot(abs(dVdRpostsample[,,"dvobs"] - gpsmooth_ss$dvobs) < 1e-10)


vdRpostcurve <- getMeanDerivCurve(x=fn.sim$time, y.mat=gpsmooth_ss$rtrue, dy.mat=gpsmooth_ss$drobs, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

vdVpostcurve <- getMeanDerivCurve(x=fn.sim$time, y.mat=gpsmooth_ss$vtrue, dy.mat=gpsmooth_ss$dvobs, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$vphi, gamma.mat=gpsmooth_ss$gamma)

id.max <- which.max(gpsmooth_ss$lp__)
# reporting and plot for true value ----

# "best" log-likelihood based on truth. phi vector found by optim with other inputs set at truth
loglik( data.matrix(fn.true[seq(1,401,length=41),c("Vtrue","Rtrue")]), init$abc, 
        c(init$vphi, init$rphi), 0.1,  fn.sim[,1:2], r)
gpsmooth_ss$lp__[id.max]
gpsmooth_ss$rphi[id.max,]
gpsmooth_ss$vphi[id.max,]
gpsmooth_ss$abc[id.max,]
gpsmooth_ss$sigma[id.max]

pdf("Gaussian Process behavior at true parameter.pdf", width = 8, height = 8)

matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1),
        ylab="R", main="at true parm")
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$rtrue), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdRpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$drobs), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdRpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
        ylab="V", main="at true parm")
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$vtrue), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdVpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$dvobs), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdVpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
dev.off()

#### Gaussian process fitting without ODE as initial value ####
gpfit <- stan(file="stan/gp-initialfit.stan",
              data=list(N=nrow(fn.sim),
                        robs=fn.sim$Rtrue,
                        vobs=fn.sim$Vtrue,
                        time=fn.sim$time),
              iter=600, chains=7, warmup = 100, cores=7)
gpfit_ss <- extract(gpfit, permuted=TRUE)

# looking at MAP - colored blue, posterior mean is colored green

hyperreta <- cbind(mean=colMeans(gpfit_ss$reta),
                   sd=apply(gpfit_ss$reta,2,sd))
hyperveta <- cbind(colMeans(gpfit_ss$veta),
                   apply(gpfit_ss$veta,2,sd))

post.noODE <- summary.post.noODE("results/STAN-noODE-noise-0.1.pdf", fn.true, fn.sim, gpfit_ss, init)

#### real simulation ####
gpsmooth1 <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time,
                           hyperparm=c(gpfit.post[,"mean"], gpfit.post[,"sd"]),
                           hyperreta=hyperreta,
                           hyperveta=hyperveta,
                           ubsigma=100),
                 iter=100, chains=5, warmup = 50, cores=5) #init = list(init,init.map,init.epost,init.marmode,"random","0",list())

init.rand <- lapply(1:10, function(dummy){
  id.max <- sample(nrow(gpfit_ss$reta),1)
  list(
    reta = gpfit_ss$reta[id.max,],
    veta = gpfit_ss$veta[id.max,],
    rphi = gpfit_ss$rphi[id.max,],
    vphi = gpfit_ss$vphi[id.max,],
    sigma = gpfit_ss$sigma[id.max]
  )
})
  
init.full <- c(post.noODE[c("init.epost","init.map","init.marmode")],
               list(init,init.rand[[1]],init.rand[[2]],init.rand[[3]],init.rand[[4]]))

gpsmooth2 <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time,
                           hyperparm=hyperparm0,
                           hyperreta=hyperreta0,
                           hyperveta=hyperveta0,
                           ubsigma=0.22),
                 init = init.full, iter=50, chains=8, warmup = 20, cores=8) #

init.pR <- init.full
init.pR <- lapply(init.pR, function(x) {x$veta <- NULL; x$vphi <- NULL; x})
gpsmooth.pR <- stan(file="stan/gp-smooth-partial-robs.stan",
                  data=list(N=nrow(fn.sim),
                            robs=fn.sim$Rtrue,
                            time=fn.sim$time,
                            hyperparm=hyperparm0,
                            hyperreta=hyperreta0,
                            hyperveta=hyperveta0,
                            ubsigma=0.22),
                  init = init.pR,
                  iter=100, chains=8, warmup = 50, cores=8) 

init.pV <- init.full
init.pV <- lapply(init.pV, function(x) {x$reta <- NULL; x$rphi <- NULL; x})
gpsmooth.pV <- stan(file="stan/gp-smooth-partial-vobs.stan",
                    data=list(N=nrow(fn.sim),
                              vobs=fn.sim$Vtrue,
                              time=fn.sim$time,
                              hyperparm=hyperparm0,
                              hyperreta=hyperreta0,
                              hyperveta=hyperveta0,
                              ubsigma=0.25),
                    init = init.pV,
                    iter=50, chains=8, warmup = 20, cores=8) #

plot.post.samples("results/STAN-ode-noise-0.1.pdf", fn.true, fn.sim, extract(gpsmooth2, permuted=TRUE), init)
plot.post.samples("results/STAN-ode-noise-0.1-partial-robs.pdf", fn.true, fn.sim, extract(gpsmooth.pR, permuted=TRUE), init)
plot.post.samples("results/STAN-ode-noise-0.1-partial-vobs.pdf", fn.true, fn.sim, extract(gpsmooth.pV, permuted=TRUE), init)

#### ad hoc analysis ####
lglik <- lapply(1:length(gpsmooth_ss$lp__), function(it){
  loglik( cbind(gpsmooth_ss$vtrue[it,], gpsmooth_ss$rtrue[it,]), 
          gpsmooth_ss$abc[it,], 
          c(gpsmooth_ss$vphi[it,],gpsmooth_ss$rphi[it,]), 
          0.1,
          fn.sim[,1:2], 
          r)
})



loglik(x = matrix(0, ncol=2,nrow=nrow(fn.sim)), 
       theta = c(0, 0.2, 1.0), 
       phi = c(1e-10, 1e+01, 1e-09, 1e+01), 
       sigma = 1.142141,
       y = fn.sim[,1:2], 
       r = r)

summary(gpsmooth_ss$lp__[c(1:50,101:150)])
# log posterior for converged chain
summary(gpsmooth_ss$lp__[-c(1:50,101:150)])
# log posterior for zero chain

# try the hard boundary on parameter as well. 