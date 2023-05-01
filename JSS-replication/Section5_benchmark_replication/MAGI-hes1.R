# Run MAGI method on Hes1 model

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) == 0){
  seed = 1
}else{
  seed = args
}

dir.create("results/magi", showWarnings = FALSE)
pdf(paste0("results/magi/Hes1-log-",seed,".pdf"))
sink(paste0("results/magi/Hes1-log-",seed,".txt"))

print(seed)
set.seed(seed)

library(magi)

source("setup-sim.R")

StartTime <- proc.time()[3]

compnames <- c("P", "M", "H")
matplot(x[, "time"], x[, -1], type = "l", lty = 1,
        xlab = "Time (min)", ylab = "Level")
matplot(y$time, y[,-1], type = "p", col = 1:(ncol(y)-1), pch = 20, add = TRUE)
legend("topright", compnames, lty = 1, col = c("black", "red", "green"))


hes1logmodel <- list(
  fOde = hes1logmodelODE,
  fOdeDx = hes1logmodelDx,
  fOdeDtheta = hes1logmodelDtheta,
  thetaLowerBound = rep(0, 7),
  thetaUpperBound = rep(Inf, 7)
)

hes1result <- MagiSolver(y, hes1logmodel, 
                         control = list(sigma = param.true$sigma, useFixedSigma = TRUE))

theta.est <- apply(hes1result$theta, 2, mean)
x0.est <- apply(hes1result$xsampled[,1,], 2, mean)

TimeUsed <- proc.time()[3] - StartTime

#### Calculate trajectory RMSE at the 33 time points
xdesolve_recon <- deSolve::ode(y = x0.est, times = seq(0, 240, by = 7.5), func = logmodelODE, parms = theta.est)

rmse_log <- sqrt(colMeans((x[match(xdesolve_recon[,1], x[,1]),-1] - xdesolve_recon[,-1])^2))
rmse_orig <- sqrt(colMeans((exp(x[match(xdesolve_recon[,1], x[,1]),-1]) - exp(xdesolve_recon[,-1]))^2))

rmse_log
rmse_orig

save.image(file = paste0("results/magi/Hes1-log-",seed,".rda"))
