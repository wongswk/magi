x <- cbind(1:4, 4:1, sin(1:4))
theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
tvec <- 1:4

dynamicalModelList <- list(
  modelODE = magi:::hes1modelODE,
  modelDtheta = magi:::hes1modelDtheta,
  modelDx = magi:::hes1modelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "Hes1 system", x, theta))


dtheta_byhand = magi:::hes1modelDtheta(theta, x, tvec)
dx_byhand = magi:::hes1modelDx(theta, x, tvec)


library("magi")
library(torch)


hes1modelODE_grad <- function(theta, x, tvec) {
  # Convert input arguments to torch tensors with requires_grad = TRUE
  theta <- torch_tensor(theta, requires_grad = TRUE)
  x <- torch_tensor(x, requires_grad = TRUE)
  tvec <- torch_tensor(tvec)
  
  # Define the original function using torch operations
  hes1modelODE_torch <- function(theta, x) {
    P = x[, 1]
    M = x[, 2]
    H = x[, 3]
    
    PMHdt = torch_empty(dim(x))
    PMHdt[, 1] = -theta[1] * P * H + theta[2] * M - theta[3] * P
    PMHdt[, 2] = -theta[4] * M + theta[5] / (1 + P^2)
    PMHdt[, 3] = -theta[1] * P * H + theta[6] / (1 + P^2) - theta[7] * H
    
    PMHdt
  }
  
  # Calculate output using torch operations
  output = hes1modelODE_torch(theta, x)
  
  # Initialize gradient matrices
  dPMHdt_dtheta = array(dim=c(nrow(output), length(theta), ncol(output)))
  dPMHdt_dx = array(dim=c(nrow(output), ncol(output), ncol(output)))
  
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
      
      output[i, j]$backward(keep_graph=TRUE)
      
      dPMHdt_dtheta[i, , j] = as_array(theta$grad)
      dPMHdt_dx[i, , j] = as_array(x$grad[i,])
    }
  }
  
  list(dPMHdt_dtheta = dPMHdt_dtheta, dPMHdt_dx = dPMHdt_dx)
}


augograd = hes1modelODE_grad(theta, x, tvec)
sum(abs(augograd$dPMHdt_dtheta - dtheta_byhand))
sum(abs(augograd$dPMHdt_dx - dx_byhand))

