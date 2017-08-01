#### utility functions ####
getMeanCurve <- function(x, y, x.new, phi.mat, delta = 1e-9, sigma.mat){
  tvec <- c(x.new,x)
  
  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2
  
  signr <- -sign(foo)
  
  t(sapply(1:nrow(phi.mat), function(it){
    sigma <- sigma.mat[it]
    phi <- phi.mat[it,]
    
    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    diag(C) <- diag(C)+delta
    diag(C)[-(1:length(x.new))] <- diag(C)[-(1:length(x.new))]+sigma
    C[1:length(x.new),-(1:length(x.new))]%*%solve(C[-(1:length(x.new)),-(1:length(x.new))], y)
  }))
}

#### start of code ####

library(rstan)
fn.true <- read.csv("data/FN.csv")
fn.true$time <- seq(0,20,0.05)
matplot(fn.true$time, data.matrix(fn.true[,-3]), type="l", lty=1)

fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=0.1)


# fit_fit <- stan(file="stan/ode-fit.stan",
#                 data=list(N=nrow(fn.sim),
#                           robs=fn.sim$Rtrue,
#                           vobs=fn.sim$Vtrue),
#                 iter=200, chains=5)

fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=41),]
matplot(fn.sim$time, data.matrix(fn.sim[,-3]), type="l", lty=1)

gpsmooth <- stan(file="stan/gp-smooth.stan",
                 data=list(N=nrow(fn.sim),
                           robs=fn.sim$Rtrue,
                           vobs=fn.sim$Vtrue,
                           time=fn.sim$time),
                 iter=100, chains=1)

traceplot(gpsmooth)
gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

plot(gpsmooth_ss$sigma, type="l",main="sigma_sq")
abline(h=0.1, col=2)
hist(gpsmooth_ss$sigma, breaks = 50,main="sigma_sq")
abline(v=0.1, col=2)

rmcurve <- getMeanCurve(x=fn.sim$time, y=fn.sim$Rtrue, x.new=fn.true$time,
                       sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$rphi)
vmcurve <- getMeanCurve(x=fn.sim$time, y=fn.sim$Vtrue, x.new=fn.true$time,
                        sigma.mat = gpsmooth_ss$sigma, phi.mat = gpsmooth_ss$vphi)

matplot(fn.sim$time, data.matrix(fn.sim[,-3]), type="l", lty=1)
matplot(fn.true$time, t(rmcurve), col="pink",add=TRUE, type="l",lty=1)
matplot(fn.true$time, t(vmcurve), col="grey",add=TRUE, type="l",lty=1)
points(fn.sim$time, fn.sim$Vtrue, col=1)
points(fn.sim$time, fn.sim$Rtrue, col=2)


