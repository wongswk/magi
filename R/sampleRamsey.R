library(deSolve)
library(CollocInfer)
source("helper/utilities.r")
source("helper/basic_hmc.R")

fhn.ode <- function(times, x, p) {
  dx <- x
  dimnames(dx) <- dimnames(x)
  dx["V"] <- p["c"] * (x["V"] - x["V"]^3 / 3 + x["R"])
  dx["R"] <- - (x["V"] - p["a"] + p["b"] * x["R"]) / p["c"]
  return(list(dx))
}
FhNvarnames <- c("V", "R")
FhNparnames <- c("a", "b", "c")

x0 <- c(-1, 1)
names(x0) <- FhNvarnames

FhNpars <- c(0.2, 0.2, 3)
names(FhNpars) <- FhNparnames

fhn.fun <- function(times, x, p, more) {
 dx <- x
 dx[, "V"] <- p["c"] * (x[, "V"] - x[, "V"]^3 / 3 + x[, "R"])
 dx[, "R"] <- -(x[, "V"] - p["a"] + p["b"] * x[, "R"]) / p["c"]
 return(dx)
}

FhNtimes <- seq(0, 20, 0.1)
FhNn <- length(FhNtimes)
out <- lsoda(x0, times = FhNtimes, fhn.ode, FhNpars)
FhNdata <- out[, 2:3] + 0.5 * matrix(rnorm(2 * FhNn), FhNn, 2)   # noisy data


# set up B-splines.
FhNrange <- c(0, 20)
breaks <- seq(0, 20, 0.1)
FhNbasis <- create.bspline.basis(range = FhNrange, norder = 4, breaks = breaks)

FhNfdPar <- fdPar(FhNbasis, int2Lfd(2), 1)
lambda <- 1000  # start pars at truth so this OK, Ramsey method's estimates converge as lambda --> infty
# initial smooth
DEfd0 <- smooth.basis(FhNtimes, FhNdata, FhNfdPar)$fd
coefs0 <- DEfd0$coef

profile.obj <- LS.setup(pars = FhNpars, fn = fhn.fun, lambda = lambda,
 times = FhNtimes, coefs = coefs0, basisvals = FhNbasis)

proc <- profile.obj$proc
lik <- profile.obj$lik

Ores2.2 <- Profile.LS(fhn.fun, FhNdata, FhNtimes, FhNpars, coefs0, FhNbasis, lambda)    # Ramsey optimization, use this as starting values for sampling

Ores2.2$pars

################################


# assume current pars (curpars) defined
Cblock <- function(q, grad=FALSE) {

    data<- FhNdata
    times <- FhNtimes
    pars <- curpars

    coefs2 = matrix(q, ncol(lik$bvals), length(q)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2, proc$bvals, pars, proc$more)

  ret <- -f
  
  if(grad) {
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals, 
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars, proc$more)
    g = as.vector(g)    
    attr(ret,"grad") <- -g
  }
  
  return(ret)
}


# assume current spline coefficients (curcoefs) are defined
parblock <- function(q, grad=FALSE) {

    data<- FhNdata
    times <- FhNtimes
    coefs <- curcoefs
    pars <- q

    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2, proc$bvals, pars, proc$more)

  ret <- -f
  
  if(grad) {
    attr(ret,"grad") <- -SplineCoefsDP(coefs, times, data, lik, proc, pars)
  }
  
  return(ret)
}


#### actual sampling

n.iter <- 10000
th.pars <- matrix(NA,n.iter,3)  
th.C <- matrix(NA,n.iter,406) 
th.pars[1,] <- Ores2.2$pars
th.C[1,] <- c(Ores2.2$coefs[,1],Ores2.2$coefs[,2])
curpars <- th.pars[1,]
curcoefs <- th.C[1,]
lliklist <- c()
lliklist[1] <- parblock(Ores2.2$pars)
accepts <- matrix(NA,n.iter,2)
accepts[1,] <- c(1,1)


for (t in 2:n.iter) {
  
  if (t %% 100 == 0) { show(c(t, lliklist[t-1], apply(accepts[1:(t-1),],2,mean))) }
  
  # Update pars
  foo <- basic_hmc(parblock, step=runif(1,0.015,0.03)* c(1,1,2.5), nsteps= 20, initial=th.pars[t-1,], return.traj = T)
  th.pars[t,] <- foo$final
  curpars <- th.pars[t,]
  accepts[t,1] <- foo$acc

  # Update C
  foo <- basic_hmc(Cblock, step=runif(1,0.005,0.01), nsteps= 20, initial=th.C[t-1,], return.traj = T)
  th.C[t,] <- foo$final
  curcoefs <- th.C[t,]
  accepts[t,2] <- foo$acc

  lliklist[t] <- foo$lpr

}


plot(lliklist)

par(mfrow=c(1,3))
hist(th.pars[,1],xlab="a") ; abline(v=0.2, col="blue", lwd=2)
hist(th.pars[,2],xlab="b") ; abline(v=0.2, col="blue", lwd=2)
hist(th.pars[,3],xlab="c") ; abline(v=3, col="blue", lwd=2)

par(mfrow=c(1,3))
plot(th.pars[,1],main="a") ; abline(h=0.2, col="blue", lwd=2)
plot(th.pars[,2],main="b") ; abline(h=0.2, col="blue", lwd=2)
plot(th.pars[,3],main="c") ; abline(h=3, col="blue", lwd=2)



#save(th.pars,th.C,lliklist,FhNdata, file="c:/data/projects/github/dynamic-systems/R/ramseySamp1.rda")