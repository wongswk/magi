# Each of the methods calls this script to set up the same simulated data

# Generate simulated dataset
hes1modelODE <- function(theta, x, tvec) {
  P = x[, 1]
  M = x[, 2]
  H = x[, 3]
  
  PMHdt = array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * P * H + theta[2] * M - theta[3] * P
  PMHdt[, 2] = -theta[4] * M + theta[5] / (1 + P^2)
  PMHdt[, 3] = -theta[1] * P * H + theta[6] / (1 + P^2) - theta[7] * H
  
  PMHdt
}

hes1logmodelODE <- function (theta, x, tvec) {
  P = exp(x[, 1])
  M = exp(x[, 2])
  H = exp(x[, 3])
  
  PMHdt <- array(0, c(nrow(x), ncol(x)))
  PMHdt[, 1] = -theta[1] * H + theta[2] * M / P - theta[3]
  PMHdt[, 2] = -theta[4] + theta[5] / (1 + P^2) / M
  PMHdt[, 3] = -theta[1] * P + theta[6] / (1 + P^2) / H - theta[7]
  
  PMHdt
}

param.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.439, 2.037, 17.904),
  sigma = c(0.15, 0.15, NA))

modelODE <- function(tvec, state, parameters) {
  list(as.vector(hes1modelODE(parameters, t(state), tvec)))
}

logmodelODE <- function(tvec, state, parameters) {
  list(as.vector(hes1logmodelODE(parameters, t(state), tvec)))
}

x <- deSolve::ode(y = param.true$x0, times = seq(0, 60 * 4, by = 0.01),
                  func = modelODE, parms = param.true$theta)

y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5), ])
names(y) <- c("time", "P", "M", "H")
y$P <- y$P * exp(rnorm(nrow(y), sd = param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd = param.true$sigma[2]))

y$H <- NaN
y$P[y$time %in% seq(7.5, 240, by = 15)] <- NaN
y$M[y$time %in% seq(0, 240, by = 15)] <- NaN

# Log transform
x[, names(y) != "time"] <- log(x[, names(y) != "time"])
y[, names(y) != "time"] <- log(y[, names(y) != "time"])
