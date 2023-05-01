# Appendix E: Demo of using autograd with MAGI

library("magi")
library("torch")
# > sessionInfo()
# other attached packages:
#   [1] torch_0.10.0 magi_1.2.2  

x <- cbind(1:4, 4:1, sin(1:4))
theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
tvec <- 1:4

dtheta_byhand = magi:::hes1logmodelDtheta(theta, x, tvec)
dx_byhand = magi:::hes1logmodelDx(theta, x, tvec)


# Define the original function using torch operations
hes1logmodelODE_torch <- function (theta, x, tvec) {
  P = exp(x[, 1])
  M = exp(x[, 2])
  H = exp(x[, 3])
  
  PMHdt <- torch_empty(dim(x))
  PMHdt[, 1] = -theta[1] * H + theta[2] * M / P - theta[3]
  PMHdt[, 2] = -theta[4] + theta[5] / (1 + P^2) / M
  PMHdt[, 3] = -theta[1] * P + theta[6] / (1 + P^2) / H - theta[7]
  
  PMHdt
}

ode_autograd <- function(ode_func_torch, theta, x, tvec) {
  # Convert input arguments to torch tensors with requires_grad = TRUE
  theta <- torch_tensor(theta, requires_grad = TRUE)
  x <- torch_tensor(x, requires_grad = TRUE)
  tvec <- torch_tensor(tvec)
  
  
  # Calculate output using torch operations
  output = ode_func_torch(theta, x)
  
  # Initialize gradient matrices
  ode_dtheta = array(dim=c(nrow(output), length(theta), ncol(output)))
  ode_dx = array(dim=c(nrow(output), ncol(output), ncol(output)))
  
  # Calculate gradients for each element in the output
  for (i in 1:nrow(output)) {
    for (j in 1:ncol(output)) {
      # Zero out gradients from previous iterations
      if (length(theta$grad) > 0) {
        theta$grad$zero_()
      }
      if (length(x$grad) > 0) {
        x$grad$zero_()
      }
      
      # use keep_graph=TRUE or retain_graph=TRUE depending on torch version
      output[i, j]$backward(retain_graph=TRUE)  
      
      ode_dtheta[i, , j] = as_array(theta$grad)
      ode_dx[i, , j] = as_array(x$grad[i,])
    }
  }
  
  list(ode_dtheta = ode_dtheta, ode_dx = ode_dx)
}

# Test for autograd correctness compared to analytic gradients
autograd = ode_autograd(hes1logmodelODE_torch, theta, x, tvec)
sum(abs(autograd$ode_dtheta - dtheta_byhand))
sum(abs(autograd$ode_dx - dx_byhand))

# Hes1 example

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

param.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(1.439, 2.037, 17.904),
  sigma = c(0.15, 0.15, NA))

modelODE <- function(tvec, state, parameters) {
  list(as.vector(hes1modelODE(parameters, t(state), tvec)))
}

x <- deSolve::ode(y = param.true$x0, times = seq(0, 60 * 4, by = 0.01),
                  func = modelODE, parms = param.true$theta)


set.seed(12321)
y <- as.data.frame(x[ x[, "time"] %in% seq(0, 240, by = 7.5), ])
names(y) <- c("time", "P", "M", "H")
y$P <- y$P * exp(rnorm(nrow(y), sd = param.true$sigma[1]))
y$M <- y$M * exp(rnorm(nrow(y), sd = param.true$sigma[2]))


y$H <- NaN
y$P[y$time %in% seq(7.5, 240, by = 15)] <- NaN
y$M[y$time %in% seq(0, 240, by = 15)] <- NaN

y.tilde <- y
y.tilde[, names(y.tilde) != "time"] <- log(y.tilde[, names(y.tilde) != "time"])

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

hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = function(theta, x, tvec) 
    ode_autograd(hes1logmodelODE_torch, theta, x, tvec)$ode_dx,
  fOdeDtheta = function(theta, x, tvec) 
    ode_autograd(hes1logmodelODE_torch, theta, x, tvec)$ode_dtheta,
  thetaLowerBound = rep(0, 7),
  thetaUpperBound = rep(Inf, 7)
)

yTest <- matrix(runif(nrow(y.tilde) * (ncol(y.tilde) - 1)),
                nrow = nrow(y.tilde), ncol = ncol(y.tilde) - 1)
thetaTest <- runif(7)
testDynamicalModel(hes1logmodel$fOde, hes1logmodel$fOdeDx, hes1logmodel$fOdeDtheta, 
                   "Hes1 log", yTest, thetaTest, y[, "time"])

# Run MagiSolver using the autograd setup (201 HMC iterations for demo only: computation is very slow)
hes1result <- MagiSolver(y.tilde, hes1logmodel, 
                         control = list(sigma = param.true$sigma, useFixedSigma = TRUE, niterHmc=201))