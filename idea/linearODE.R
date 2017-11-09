library(deSolve)

# set up linear ODE --------------------------------------------------------

parameters <- c(a=1, b=-0.2)
state <- c(X = 1)

linearModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Xdt <- a + b * X
    list(c(Xdt))
  })
}

times <- seq(0, 20, by = 0.05)
out <- ode(y = state, times = times, func = linearModel, parms = parameters)
plot(out[,"time"], out[, "X"], type="l")


# analytical solution for ODE ---------------------------------------------


XfunAnalytic <- function(t, par){
  par[1] + par[2] * exp(par[3] * t)
}

par.true <- with(as.list(c(state, parameters)), {
  c(-a/b, X+a/b, b)
})


Xfunc <- function(t) XfunAnalytic(t, par.true)
plot.function(Xfunc, 0, 20, col=2, add=TRUE)


# simulate observations --------------------------------------------------
nobs <- 21
tI <- seq(0, 10, length=nobs)
y <- rnorm(length(tI))*0.1 + Xfunc(tI)

points(tI, y)

# frequentist inference ------------------------------------------------------

mse <- function(par){
  mean((y - XfunAnalytic(tI, par))^2)
}

frequentist <- optim(par.true, mse, method = "BFGS")
sigmaEst <- sqrt(frequentist$value)

llik <- function(par){
  sum(dnorm(y, XfunAnalytic(tI, par), sigmaEst, log = TRUE))
}

frequentist <- optim(frequentist$par, llik, method = "BFGS", hessian = TRUE, 
                     control = list(fnscale = -1))

parEst <- frequentist$par
varMat <- solve(-frequentist$hessian)
sqrt(diag(varMat))

(parEst - par.true)/sqrt(diag(varMat))
