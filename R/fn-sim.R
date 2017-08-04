#### start of code ####
source("R/func.R")
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
gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time),
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
              iter=200, chains=7, warmup = 100, cores=7)
gpfit_ss <- extract(gpfit, permuted=TRUE)

id.plot <- seq(1,length(gpfit_ss$lp__),length=100)
id.plot <- unique(as.integer(id.plot))

# looking at MAP - colored blue, posterior mean is colored green
id.max <- which.max(gpfit_ss$lp__)
# looking at posterior distribution

pdf("GP fit (no ODE information).pdf", width = 8, height = 8)
layout(1)
matplot(fn.true$time, data.matrix(fn.true[,c(1:2)]), type="l", lty=1, col=c(2,1), 
        ylab="R & V", main="full posterior")
points(fn.sim$time, fn.sim$Rtrue, col=1)
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpfit_ss$rtrue[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.sim$time, t(gpfit_ss$vtrue[id.plot,]), col="pink",add=TRUE, type="p",lty=1, pch=20)

matplot(fn.true$time, data.matrix(fn.true[,c(1:2)]), type="l", lty=1, col=c(2,1), 
        ylab="R & V", main="full posterior", add=TRUE)

lines(fn.sim$time, colMeans(gpfit_ss$rtrue), col=3)
lines(fn.sim$time, colMeans(gpfit_ss$vtrue), col=3)

lines(fn.sim$time, gpfit_ss$rtrue[id.max,], col=4)
lines(fn.sim$time, gpfit_ss$vtrue[id.max,], col=4)

layout(matrix(1:4,2,byrow = TRUE))
hist(gpfit_ss$rphi[,1])
abline(v=init$rphi[1], col=2)
abline(v=mean(gpfit_ss$rphi[,1]), col=3)
abline(v=gpfit_ss$rphi[id.max,1], col=4)
hist(gpfit_ss$rphi[,2])
abline(v=init$rphi[2], col=2)
abline(v=mean(gpfit_ss$rphi[,2]), col=3)
abline(v=gpfit_ss$rphi[id.max,2], col=4)
hist(gpfit_ss$vphi[,1])
abline(v=init$vphi[1], col=2)
abline(v=mean(gpfit_ss$vphi[,1]), col=3)
abline(v=gpfit_ss$vphi[id.max,1], col=4)
hist(gpfit_ss$vphi[,2])
abline(v=init$vphi[2], col=2)
abline(v=mean(gpfit_ss$vphi[,2]), col=3)
abline(v=gpfit_ss$vphi[id.max,2], col=4)

layout(1)
hist(gpfit_ss$sigma)
abline(v=init$sigma, col=2)
abline(v=mean(gpfit_ss$sigma), col=3)
abline(v=gpfit_ss$sigma[id.max], col=4)
dev.off()

init.map <- list(
  reta = gpfit_ss$reta[id.max,],
  veta = gpfit_ss$veta[id.max,],
  rphi = gpfit_ss$rphi[id.max,],
  vphi = gpfit_ss$vphi[id.max,],
  sigma = gpfit_ss$sigma[id.max]
)

init.epost <- list(
  reta = colMeans(gpfit_ss$reta),
  veta = colMeans(gpfit_ss$veta),
  rphi = colMeans(gpfit_ss$rphi),
  vphi = colMeans(gpfit_ss$vphi),
  sigma = mean(gpfit_ss$sigma)
)

#### real simulation ####
gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time),
                 iter=100, chains=3, warmup = 50, cores=3, 
                 init = list(init,init.map,init.epost))


gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)
max(gpsmooth_ss$lp__)
id.max <- which.max(gpsmooth_ss$lp__)

plot(gpsmooth_ss$sigma, type="l",main="sigma")
abline(h=0.1, col=2)
hist(gpsmooth_ss$sigma, breaks = 50,main="sigma")
abline(v=0.1, col=2)

vdRpostcurve <- getMeanDerivCurve(x=fn.sim$time, y.mat=gpsmooth_ss$rtrue, dy.mat=gpsmooth_ss$drobs, x.new=fn.true$time,
                                  sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi, gamma.mat=gpsmooth_ss$gamma)

vdVpostcurve <- getMeanDerivCurve(x=fn.sim$time, y.mat=gpsmooth_ss$vtrue, dy.mat=gpsmooth_ss$dvobs, x.new=fn.true$time,
                                  sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$vphi, gamma.mat=gpsmooth_ss$gamma)

lglik <- lapply(1:length(gpsmooth_ss$lp__), function(it){
  loglik( cbind(gpsmooth_ss$vtrue[it,], gpsmooth_ss$rtrue[it,]), 
          gpsmooth_ss$abc[it,], 
          c(gpsmooth_ss$vphi[it,],gpsmooth_ss$rphi[it,]), 
          0.1,
          fn.sim[,1:2], 
          r)
})

#### reporting and plot for real simulation ####

summary(unlist(lglik))
summary(gpsmooth_ss$lp__)
lglik[[id.max]]
gpsmooth_ss$rphi[id.max,]
gpsmooth_ss$vphi[id.max,]
gpsmooth_ss$abc[id.max,]
gpsmooth_ss$sigma[id.max]


pdf("MCMC initialize different values (vaguely informative prior).pdf", width = 8, height = 8)
id.plot <- seq(1,nrow(gpsmooth_ss$abc),length=100)
id.plot <- unique(as.integer(id.plot))

traceplot(gpsmooth, pars=names(gpsmooth)[1:8])

matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1), 
        ylab="R", main="full posterior")
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$rtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdRpostcurve[id.plot,]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$drobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdRpostcurve[id.plot,]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
        ylab="V", main="full posterior")
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$vtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdVpostcurve[id.plot,]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$dvobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdVpostcurve[id.plot,]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)




layout(matrix(1:4,2,byrow = TRUE))
hist(gpsmooth_ss$abc[,1], main="a")
abline(v=init$abc[1], col=2)
hist(gpsmooth_ss$abc[,2], main="b")
abline(v=init$abc[2], col=2)
hist(gpsmooth_ss$abc[,3], main="c")
abline(v=init$abc[3], col=2)

hist(gpsmooth_ss$sigma, main="sigma")
abline(v=init$sigma, col=2)

hist(gpsmooth_ss$rphi[,1], main="phi1_R")
hist(gpsmooth_ss$rphi[,2], main="phi2_R")
hist(gpsmooth_ss$vphi[,1], main="phi1_V")
hist(gpsmooth_ss$vphi[,2], main="phi2_V")

layout(1)

#' discrepency between STAN log posterior and log likelihood
#' need to implement the model in plain R later
plot(unlist(lglik), gpsmooth_ss$lp__, xlab="log likelihood from Rscript",
     ylab="log posterior from STAN")


matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1),
        ylab="R", main="maximum a posterior")
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$rtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdRpostcurve[id.max,,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$drobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdRpostcurve[id.max,,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)


matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
        ylab="V", main="maximum a posterior")
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpsmooth_ss$vtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdVpostcurve[id.max,,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpsmooth_ss$dvobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdVpostcurve[id.max,,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
dev.off()

save(gpsmooth_ss, vdRpostcurve, vdVpostcurve, file="dump.RData")

loglik(x = matrix(0, ncol=2,nrow=nrow(fn.sim)), 
       theta = c(0, 0.2, 1.0), 
       phi = c(1e-10, 1e+01, 1e-09, 1e+01), 
       sigma = 1.142141,
       y = fn.sim[,1:2], 
       r = r)
