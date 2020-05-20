library(gpds)

# remove results that don't have common seed
rdaDir <- "/Volumes/TimeMachineBackup/Workspace/DynamicSys/results/cpp/7param-917321ad09779c5bd6932651b8291acc4c236d22/"
subdirs <- list.dirs(rdaDir)[-1]

env_all <- list()
env_all$notemper <- new.env()
env_all$warmstart <- new.env()
env_all$updatephi <- new.env()
# env_all$reoptimizephi <- new.env()

load(paste0(rdaDir, "/variablephi-notemper/hes1log_summary.rda"), envir = env_all$notemper)
load(paste0(rdaDir, "/variablephi-temper-warmstart/hes1log_summary.rda"), envir = env_all$warmstart)
load(paste0(rdaDir, "/variablephi-temper-warmstart-updatephi/hes1log_summary.rda"), envir = env_all$updatephi)
# load(paste0(rdaDir, "/variablephi-temper-warmstart-reoptimizephi/hes1log_summary.rda"), envir = env_all$reoptimizephi)

all_seeds <- sapply(env_all, function(x) gsub(".*log-([0-9]+)-7param.*", "\\1", x$rda_files))
for (i in 2:length(env_all))
  stopifnot(all(all_seeds[,1] == all_seeds[,i]))

rmse.obs.all <- sapply(env_all$notemper$ours, function(x) x$rmseOdePM)
round(rowMeans(rmse.obs.all), digits=4)
hist(rmse.obs.all[1,], breaks=100, main="RMSE on observation without tempering")
abline(v=0.3)
mtext("component 1: P")
hist(rmse.obs.all[2,], breaks=100, main="RMSE on observation without tempering")
abline(v=0.3)
mtext("component 2: M")

id_to_update <- which((rmse.obs.all[1,] > 0.3) | (rmse.obs.all[2,] > 0.3))

ours <- env_all$warmstart$ours
ours[id_to_update] <- env_all$updatephi$ours[id_to_update]

oursPostExpX <- env_all$warmstart$oursPostExpX
oursPostExpX[id_to_update] <- env_all$updatephi$oursPostExpX[id_to_update]

oursExpXdesolvePM <- env_all$warmstart$oursExpXdesolvePM
oursExpXdesolvePM[id_to_update] <- env_all$updatephi$oursExpXdesolvePM[id_to_update]


load(paste0(rdaDir, "/variablephi-notemper/", env_all$notemper$rda_files[1]), envir = .GlobalEnv)
rm(env_all)
save.image(paste0(rdaDir, "Hes1Hybri.rda"))

rowId <- sapply(xsim$time, function(x) which(abs(x-times) < 1e-6))
for (i in 1:length(ours)) {
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  ours[[i]]$rmseOdeExpPM <- sqrt(apply((exp(xdesolvePM.obs) - exp(xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
rmse_orig <- round(apply(sapply(ours, function(x) x$rmseOdeExpPM), 1, mean), digits=4)
print(rmse_orig)

oursPostExpX <- sapply(oursPostExpX, identity, simplify = "array")
oursExpXdesolvePM <- sapply(oursExpXdesolvePM, function(x) x[,-1], simplify = "array")
xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])

ylim_lower <- c(1.5, 0.5, 0)
ylim_upper <- c(9.0, 3.1, 19)

pdf(width = 20, height = 5, file=paste0(rdaDir, "/posteriorExpxHes1HybridOursIllustration.pdf"))
par(mfrow=c(1, ncol(xsim)+1))

matplot(xtrue[, "time"], exp(xtrue[, -1]), type="l", lty=1, col=c(4,6,"goldenrod1"), xlab="time", ylab=NA)
matplot(xsim.obs$time, exp(xsim.obs[,-1]), type="p", col=c(4,6,"goldenrod1"), pch=20, add = TRUE)
mtext('sample observations', cex=1.5)
legend("topright", c("true P", "true M", "true H", "observed P", "observed M"), 
       lty=c(1,1,1,NA,NA), pch=c(NA,NA,NA,20,20), col=c(4,6,"goldenrod1"), cex=1.5)

phiVisualization <- rbind(
  c(2.07, 0.38, 0.45),
  c(64, 40, 23)
)
compnames <- c("P", "M", "H")

# smooth visualization with illustration
xdesolveTRUE <-ours[[1]]$xdesolveTRUE
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
id <- seq(1, nrow(xdesolveTRUE), by=50)
xdesolveTRUE <- xdesolveTRUE[id,]

ourExpXdesolveLB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.025))
ourExpXdesolveMed <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.5))
ourExpXdesolveUB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.975))


for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.975)
  
  ourEst <- exp(getMeanCurve(xsim$time, log(ourEst), xdesolveTRUE[,1], 
                             t(phiVisualization[,i]), 0, 
                             kerneltype=config$kernel, deriv = FALSE))
  ourUB <- exp(getMeanCurve(xsim$time, log(ourUB), xdesolveTRUE[,1], 
                            t(phiVisualization[,i]), 0, 
                            kerneltype=config$kernel, deriv = FALSE))
  ourLB <- exp(getMeanCurve(xsim$time, log(ourLB), xdesolveTRUE[,1], 
                            t(phiVisualization[,i]), 0, 
                            kerneltype=config$kernel, deriv = FALSE))
  
  
  times <- xdesolveTRUE[,1]
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i], cex=1.5)
  
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = "skyblue", lty = 1, density = 10, angle = -45)
  
  polygon(c(times, rev(times)), c(ourExpXdesolveUB[,i], rev(ourExpXdesolveLB[,i])),
          col = "grey80", border = "grey80", lty = 1, density = 10, angle = 45)
  
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourExpXdesolveMed[,i], lwd=3)
  lines(times, ourEst, col="forestgreen", lwd=3)
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth", "median posterior mean", "median reconstructed trajectory", 
                   "CI on posterior mean", "CI on reconstructed trajectory"), lty=c(1,1,1,0,0), lwd=c(4,3,3,0,0),
       col = c("red", "forestgreen", "black", NA, NA), density=c(NA, NA, NA, 40, 40), fill=c(0, 0, 0, "skyblue", "grey80"),
       border=c(0, 0, 0, "skyblue", "grey80"), angle=c(NA,NA,NA,-45,45), x.intersp=c(2.5,2.5,2.5,0, 0),  bty = "n", cex=1.8)
dev.off()
