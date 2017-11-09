library(deSolve)
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

XfunAnalytic <- function(t, par){
  par[1] + par[2] * exp(par[3] * t)
}

par.true <- with(as.list(c(state, parameters)), {
  c(-a/b, X+a/b, b)
})

linearModelAnalytical <- function(state, parameters) {
  par <- with(as.list(c(state, parameters)), {
    c(-a/b, X+a/b, b)
  })
  function(t){
    par[1] + par[2] * exp(par[3] * t)
  }
}

Xfunc <- linearModelAnalytical(state, parameters)
plot.function(Xfunc, 0, 20, col=2, add=TRUE)

tI <- seq(0, 10, 0.5)
y <- rnorm(length(tI))*0.1 + Xfunc(tI)

plot(tI, y)

mse <- function(par){
  mean((y - XfunAnalytic(tI, par))^2)
}

frequentist <- optim(par.true, mse, method = "BFGS", hessian = TRUE)
sigmaSqEst <- frequentist$value

