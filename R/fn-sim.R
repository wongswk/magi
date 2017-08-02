#### start of code ####
source("R/func.R")
library(rstan)
fn.true <- read.csv("data/FN.csv")
fn.true$time <- seq(0,20,0.05)
matplot(fn.true$time, data.matrix(fn.true[,-3]), type="l", lty=1)

abc = c(0.2, 0.2, 3)

fn.true$dVtrue = with(fn.true, abc[3] * (Vtrue - Vtrue^3/3.0 + Rtrue))
fn.true$dRtrue = with(fn.true, -1.0/abc[3] * (Vtrue - abc[1] + abc[2]*Rtrue))

fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=0.1)


# fit_fit <- stan(file="stan/ode-fit.stan",
#                 data=list(N=nrow(fn.sim),
#                           robs=fn.sim$Rtrue,
#                           vobs=fn.sim$Vtrue),
#                 iter=200, chains=5)
set.seed(123)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=41),]
matplot(fn.sim$time, data.matrix(fn.sim[,-3]), type="l", lty=1)

init <- list(
  abc=c(0.2,0.2,3),
  rphi=c(302.6, 2.02),
  vphi=c(145.6, 1.21)
)

tvec41 <- fn.sim$time
foo <- outer(tvec41, t(tvec41),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)


init$reta <- solve(t(chol(calCov(init$rphi, r)$C)), fn.true[seq(1,401,length=41),c("Rtrue")])
init$veta <- solve(t(chol(calCov(init$vphi, r)$C)), fn.true[seq(1,401,length=41),c("Vtrue")])

stopifnot(abs(getX(r, t(init$rphi), t(init$reta)) - fn.true[seq(1,401,length=41),c("Rtrue")]) < 1e-14)

gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time),
                 iter=2, chains=1, init=list(init), warmup = 0)


gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)
gpsmooth_ss$lp__

# init has lp value -164.7908 
# simulation around -180

stopifnot(abs(gpsmooth_ss$vtrue[1,] - fn.true[seq(1,401,length=41),c("Vtrue")])<1e-10)
stopifnot(abs(gpsmooth_ss$rtrue[1,] - fn.true[seq(1,401,length=41),c("Rtrue")])<1e-10)

covR.init <- calCov(init$rphi, r)
covV.init <- calCov(init$vphi, r)

stopifnot(abs(gpsmooth_ss$C_rphi[1,,] - covR.init$C)<1e-12)
stopifnot(abs(gpsmooth_ss$L_C_rphi[1,,] - t(chol(covR.init$C))) < 1e-12)
stopifnot(abs(gpsmooth_ss$K_rphi[1,,] - covR.init$Kphi) < 1e-10)
stopifnot(abs(gpsmooth_ss$dC_rphi[1,,] - covR.init$Cprime) < 1e-12)
stopifnot(abs(gpsmooth_ss$ddC_rphi[1,,] - covR.init$Cdoubleprime) < 1e-12)
stopifnot(abs(gpsmooth_ss$m_rphi_rtrue[1,] - covR.init$mphi%*%fn.true[seq(1,401,length=41),c("Rtrue")]) < 1e-11)
stopifnot(abs(gpsmooth_ss$drobs - fn.true[seq(1,401,length=41),c("dRtrue")])<1e-10)

gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time),
                 iter=2, chains=1, init=list(init), warmup = 0)

traceplot(gpsmooth)
gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)
max(gpsmooth_ss$lp__)
id.max <- which.max(gpsmooth_ss$lp__)

plot(gpsmooth_ss$sigma, type="l",main="sigma")
abline(h=0.1, col=2)
hist(gpsmooth_ss$sigma, breaks = 50,main="sigma")
abline(v=0.1, col=2)

plot(gpsmooth_ss$gamma, type="l",main="gamma")
hist(gpsmooth_ss$gamma, breaks = 50,main="gamma")


vdRmcurve <- getMeanDerivCurve(x=fn.sim$time, y=fn.sim$Rtrue, dy=fn.sim$dRtrue, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

vdVmcurve <- getMeanDerivCurve(x=fn.sim$time, y=fn.sim$Vtrue, dy=fn.sim$dVtrue, x.new=fn.true$time,
                               sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

Rpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$rphi, eta.mat = gpsmooth_ss$reta)
Vpostsample <- getX(r=as.matrix(dist(fn.sim$time)), phi.mat = gpsmooth_ss$vphi, eta.mat = gpsmooth_ss$veta)

dVdRpostsample <- getdVdR(abc.mat = gpsmooth_ss$abc, rtrue.mat = gpsmooth_ss$rtrue, vtrue.mat = gpsmooth_ss$vtrue)

matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(Rpostsample), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$rtrue), col="skyblue",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(dVdRpostsample[,,"drobs"]), col="skyblue",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(Vpostsample), col="skyblue",add=TRUE, type="l",lty=1)

matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1))
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(dVdRpostsample[,,"dvobs"]), col="skyblue",add=TRUE, type="l",lty=1)

save(gpsmooth_ss, vdRmcurve, vdVmcurve, file="dump.RData")

# "best" log-likelihood based on truth. phi vector found by optim with other inputs set at truth
loglik( data.matrix(fn.true[seq(1,401,length=41),c("Vtrue","Rtrue")]), init$abc, 
        c(init$vphi, init$rphi), 0.1,  fn.sim[,1:2], r)
# around 385.3999

lglik <- lapply(1:length(gpsmooth_ss$lp__), function(it){
  loglik( cbind(gpsmooth_ss$vtrue[it,], gpsmooth_ss$rtrue[it,]), 
          gpsmooth_ss$abc[it,], 
          c(gpsmooth_ss$vphi[it,],gpsmooth_ss$rphi[it,]), 
          0.1,
          fn.sim[,1:2], 
          r)
})

lglik[[id.max]]
gpsmooth_ss$abc[id.max,]
plot(unlist(lglik), gpsmooth_ss$lp__)
# around 81.03348

#' discrepency between STAN log posterior and log likelihood
#' need to implement the model in plain R later
gpsmooth_ss$rphi[id.max,]
gpsmooth_ss$vphi[id.max,]
