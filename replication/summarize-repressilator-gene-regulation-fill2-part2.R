
# sink(paste0(rdaDir,"result-partobs-fill",config$filllevel, ".txt"))

load("results/repressilator-gene-regulation-log/small-summary-partobs-fill2.rda")
load("results/repressilator-gene-regulation-log/repressilator-gene-regulation-log-1-partobs-fill2.rda")

# load(paste0(rdaDir, rda_files[1]), envir = .GlobalEnv)
rdaDir <- rdaDirSummary

library(magi)
library(xtable)

# Average the posterior mean RMSEs for the different seeds
starttime <- xsim$time[which(is.finite(xsim[,2]))[1]]
xsim.obs <- xsim.obs[xsim.obs$time >= starttime,]
rowId <- sapply(xsim.obs$time, function(x) which(abs(x-times) < 1e-6))

for (i in 1:length(ours)) {
  rowId <- sapply(xsim.obs$time, function(x) which(abs(x-ours[[i]]$xdesolveTRUE[,1]) < 1e-6))
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  rowId <- sapply(xsim.obs$time, function(x) which(abs(x-ours[[i]]$xdesolvePM[,1]) < 1e-6))
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xsim$time) < 1e-6))
  xpostmean.obs <- ours[[i]]$xpostmean[rowId,]
  
  ours[[i]]$rmseOdeExpPM <- sqrt(apply((exp(xdesolvePM.obs) - exp(xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
  ours[[i]]$rmseInferredTrajectory <- sqrt(apply((exp(xpostmean.obs) - exp(xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
rmse_reconstructed_orig <- round(apply(sapply(ours, function(x) x$rmseOdeExpPM), 1, mean), digits=4)
print(rmse_reconstructed_orig)
rmse_inferred_orig <- round(apply(sapply(ours, function(x) x$rmseInferredTrajectory), 1, mean), digits=4)
print(rmse_inferred_orig)


oursPostExpX <- sapply(oursPostExpX, identity, simplify = "array")
oursPostX <- sapply(oursPostX, identity, simplify = "array")
oursExpXdesolvePM_time <- log(oursExpXdesolvePM[[1]][,1])
oursExpXdesolvePM <- sapply(oursExpXdesolvePM, function(x) x[,-1], simplify = "array")
xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = times, func = odemodel$modelODE, parms = pram.true$theta)
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])

ylim_lower <- c(0, 0, 0, 0, 0, 0)
ylim_upper <- rep(90, 6)

pdf(width = 20, height = 5, file=paste0(rdaDir, "/PmExpX-NoNumSolver.pdf"))
# layout(rbind(c(1,2,3,4), c(5,5,5,5)), heights = c(5,1))

matplot(xtrue[, "time"], exp(xtrue[, -1]), type="l", lty=1, col=c(1:6), xlab="time", ylab=NA)
matplot(xsim.obs$time, exp(xsim.obs[,-c(1, 5:7)]), type="p", col=c(1:6), pch=20, add = TRUE)
mtext('sample observations', cex=1.5)
legend("topleft", c("true mRNA lacI", "true mRNA tetR", "true mRNA cI", "true protein lacI", "true protein tetR", "true protein cI",
                    "observed mRNA lacI", "observed mRNA tetR", "observed mRNA cI"), 
       lty=c(1,1,1,1,1,1,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,NA,20,20,20), col=c(1:6, 1:3), cex=1.5)

phiVisualization <- phiExogenous
compnames <- c("mRNA lacI", "mRNA tetR", "mRNA cI", "protein lacI", "protein tetR", "protein cI")

# smooth visualization with illustration
xdesolveTRUE <-ours[[1]]$xdesolveTRUE
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
id <- seq(1, nrow(xdesolveTRUE), by=1)
xdesolveTRUE <- xdesolveTRUE[id,]


for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.975)
  
  ourEst <- exp(magi:::getMeanCurve(xsim$time, log(ourEst), xdesolveTRUE[,1],
                                    t(phiVisualization[,i]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE))
  ourUB <- exp(magi:::getMeanCurve(xsim$time, log(ourUB), xdesolveTRUE[,1],
                                   t(phiVisualization[,i]), 0, 
                                   kerneltype=config$kernel, deriv = FALSE))
  ourLB <- exp(magi:::getMeanCurve(xsim$time, log(ourLB), xdesolveTRUE[,1],
                                   t(phiVisualization[,i]), 0, 
                                   kerneltype=config$kernel, deriv = FALSE))
  
  
  times <- xdesolveTRUE[,1]
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  if (i >= 4){
    mtext(paste(compnames[i], "component (Unobserved)"), cex=1.5)  
  }else{
    mtext(paste(compnames[i], "component (Partially Observed)"), cex=1.5)  
  }
  
  
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  
  lines(times[-1], xdesolveTRUE[-1,1+i], col="red", lwd=4)
  lines(times, ourEst, col="forestgreen", lwd=3)
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)

legend("center", c("truth", "median of all \n inferred trajectories", "95% interval from the \n 2.5 and 97.5 percentile n of all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "forestgreen", NA), fill=c(0, 0,"skyblue"), text.width=c(0, 0.4, 0.05), bty = "n",
       border=c(0, 0, "skyblue"), pch=c(NA, NA, 15), cex=1.8, horiz=FALSE)
mtext("posterior mean of Exp(X)")
dev.off()

pdf(width = 20, height = 5, file=paste0(rdaDir, "/ExpPmX.pdf"))


# smooth visualization with illustration

id <- seq(1, nrow(oursExpXdesolvePM), by=2)
ourExpXdesolveLB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.025))
ourExpXdesolveMed <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.5))
ourExpXdesolveUB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.975))
timesLUB <- oursExpXdesolvePM_time[id]

for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
  
  ourEst <- exp(magi:::getMeanCurve(xsim$time, (ourEst), xdesolveTRUE[,1],
                                    t(phiVisualization[,i]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE))
  ourUB <- exp(magi:::getMeanCurve(xsim$time, (ourUB), xdesolveTRUE[,1],
                                   t(phiVisualization[,i]), 0, 
                                   kerneltype=config$kernel, deriv = FALSE))
  ourLB <- exp(magi:::getMeanCurve(xsim$time, (ourLB), xdesolveTRUE[,1],
                                   t(phiVisualization[,i]), 0, 
                                   kerneltype=config$kernel, deriv = FALSE))
  
  
  times <- xdesolveTRUE[,1]
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = "skyblue", lty = 1, density = 10, angle = -45)
  
  polygon(c(timesLUB, rev(timesLUB)), c(ourExpXdesolveUB[,i], rev(ourExpXdesolveLB[,i])),
          col = "grey80", border = "grey80", lty = 1, density = 10, angle = 45)
  
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(timesLUB, ourExpXdesolveMed[,i], lwd=3)
  lines(times, ourEst, col="forestgreen", lwd=3)
  
  if (i >= 4){
    mtext(paste(compnames[i], "component (Unobserved)"), cex=1.5)  
  }else{
    mtext(paste(compnames[i], "component (Partially Observed)"), cex=1.5)  
  }
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = TRUE)
legend("center", c("truth", "median inferred trajectory", "median reconstructed trajectory", 
                   "CI on inferred trajectory", "CI on reconstructed trajectory"), lty=c(1,1,1,0,0), lwd=c(4,3,3,0,0),
       col = c("red", "forestgreen", "black", NA, NA), density=c(NA, NA, NA, 40, 40), fill=c(0, 0, 0, "skyblue", "grey80"), 
       border=c(0, 0, 0, "skyblue", "grey80"), angle=c(NA,NA,NA,-45,45), x.intersp=c(2.5,2.5,2.5,0, 0),  bty = "n", cex=1.8)
mtext("Exp(posterior mean of X)")
dev.off()


oursPostTheta <- sapply(oursPostTheta, identity, simplify = "array")

printr <- function(x) format(round(x, 4), nsmall=4)
tablizeEstErr <- function(est, err){
  paste(format(round(est, 4), nsmall=4), "\\pm", format(round(err, 4), nsmall=4))
}

mean_est <- rbind(
  rowMeans(oursPostTheta[,1,])
)

sd_est <- rbind(
  apply(oursPostTheta[,1,], 1, sd)
)

tab <- rbind(
  c("Ours", tablizeEstErr(mean_est[1,],sd_est[1,]))
)
tab <- data.frame(tab)
colnames(tab) <- c("Method", c("alpha0", "alpha", "n", "beta"))
rownames(tab) <- NULL
coverage <- rbind(
  printr(rowMeans((oursPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= oursPostTheta[,3,])))
)
coverage <- cbind(c("Ours"), coverage)
colnames(coverage) <- c("Method", c("alpha0", "alpha", "n", "beta"))

tab <- cbind(c("truth", pram.true$theta), t(tab), t(coverage))

theta_rmse <- sqrt(rowMeans((oursPostTheta[,1,] - pram.true$theta)^2))
tab <- cbind(tab, c("theta rmse", printr(theta_rmse)))

print("length(ours)=")
print(length(ours))
print("rmse_inferred_orig")
print(rmse_inferred_orig)
print("rmse_reconstructed_orig")
print(rmse_reconstructed_orig)
print(tab)
print(xtable(tab))
