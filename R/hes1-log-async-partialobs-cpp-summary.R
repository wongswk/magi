# Summarize the results
library(gpds)

pdf_files <- list.files("~/Workspace/DynamicSys/results/cpp/good/")
pdf_files <- pdf_files[grep("Hes1-log-.*\\.pdf", pdf_files)]
rda_files <- gsub("\\.pdf", ".rda", pdf_files)
write.table(rda_files, file="good_hes1_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
for (f in files){
  paste0(strsplit(f, "\\.")[[1]][1], ".rda")
}

rda_files <- read.table("~/Workspace/DynamicSys/good_hes1_list.txt", sep="\n")

config <- list()
config$modelName <- "Hes1-log"
config$noise <- c(0.15,0.15,0.1)
rdaDir <- "~/Workspace/DynamicSys/results/cpp/"   ## where ours rda saved

hes1logmodel <- list(
  name="Hes1-log",
  fOde=gpds:::hes1logmodelODE,
  fOdeDx=gpds:::hes1logmodelDx,
  fOdeDtheta=gpds:::hes1logmodelDtheta,
  thetaLowerBound=rep(0,7),
  thetaUpperBound=rep(Inf,7)
)

hes1modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1modelODE(parameters, t(state))))
}

## Helper function adapted from Visualization to extract trajectories and RMSE
rmsePostSamples <- function(xtrue, dotxtrue, xsim, gpode, param, config, odemodel=NULL){
  xpostmean <- apply(gpode$xsampled, 2:3, mean)
  if(!is.null(odemodel)){
    times <- sort(unique(round(c(odemodel$times, xsim$time, xtrue[,"time"]), 7)))
    xdesolveTRUE <- deSolve::ode(y = param$x0, times = times, func = odemodel$modelODE, parms = param$theta)
    
    mapId <- which.max(gpode$lglik)
    ttheta <- gpode$theta[mapId,]
    tx0 <- gpode$xsampled[mapId,1,]
    xdesolveMAP <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    ttheta <- colMeans(gpode$theta)
    tx0 <- colMeans(gpode$xsampled[,1,])
    xdesolvePM <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    rowId <- sapply(xsim$time, function(x) which(abs(x-times) < 1e-6))
    xdesolveTRUE.obs <- xdesolveTRUE[rowId,-1]
    xdesolveMAP.obs <- xdesolveMAP[rowId,-1]
    xdesolvePM.obs <- xdesolvePM[rowId,-1]
    
    xdesolveSamples <- parallel::mclapply(1:16, function(dummy){
      mapId <- sample(1:length(gpode$lglik), 1)
      ttheta <- gpode$theta[mapId,]
      tx0 <- gpode$xsampled[mapId,1,]
      deSolve::ode(y = tx0, times = odemodel$times, func = odemodel$modelODE, parms = ttheta)
    }, mc.cores = 8)
    
    rmseTrue <- sqrt(apply((xdesolveTRUE.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseWholeGpode <- sqrt(apply((xpostmean - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseOdeMAP <- sqrt(apply((xdesolveMAP.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseOdePM <- sqrt(apply((xdesolvePM.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    
    rmselist <- list(
      true = paste0(round(rmseTrue, 3), collapse = "; "),
      wholeGpode = paste0(round(rmseWholeGpode, 3), collapse = "; "),
      odeMAP = paste0(round(rmseOdeMAP, 3), collapse = "; "),
      odePM = paste0(round(rmseOdePM, 3), collapse = "; ")
    )
    
    config$rmse <- rmselist
    return( list(xdesolveTRUE = xdesolveTRUE, xdesolvePM = xdesolvePM, rmseOdePM = rmseOdePM))
  }
  
}

#### Grab the data for each seed and save in a list
ours <- list()
oursPostX <- list()
oursPostTheta <- list()
oursPostExpX <- list()
oursExpXdesolvePM <- list()

rda_files <- as.character(unlist(rda_files))
for (f in rda_files) {
  show(f)
  load(paste0(rdaDir, f))
  ours[[f]] <- rmsePostSamples(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
  oursPostX[[f]] <- cbind(
    apply(gpode$xsampled, 2:3, mean), 
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  oursPostTheta[[f]] <- cbind(
    apply(gpode$theta, 2, mean), 
    apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  )
  
  xsampledexp <- exp(gpode$xsampled)
  oursPostExpX[[f]] <- cbind(
    apply(xsampledexp, 2:3, mean), 
    apply(xsampledexp, 2:3, function(x) quantile(x, 0.025)),
    apply(xsampledexp, 2:3, function(x) quantile(x, 0.975))
  )
  
  ttheta <- colMeans(gpode$theta)
  exptx0 <- colMeans(xsampledexp[,1,])
  xdesolvePM <- deSolve::ode(y = exptx0, times = times, func = hes1modelODE, parms = ttheta)
  oursExpXdesolvePM[[f]] <- xdesolvePM
}

ramsayPostX0 <- list()
ramsayPostTheta <- list()
ramsayXdesolvePM <- list()
ramsayRmseLog <- list()
ramsayRmseOrig <- list()
for (f in rda_files){
  env_ramsay <- new.env()
  load(paste0(rdaDir, "ramsay-optimation-", f), envir = env_ramsay)
  ramsayXdesolvePM[[f]] <- env_ramsay$xdesolveRamsay
  ramsayRmseLog[[f]] <- env_ramsay$rmse_log
  ramsayRmseOrig[[f]] <- env_ramsay$rmse_orig
  ramsayPostTheta[[f]] <- env_ramsay$best.pars
  ramsayPostX0[[f]] <- env_ramsay$best.pars.x0
}
  
save.image(file=paste0(rdaDir,"hes1log_summary.rda"))

library(gpds)
library(xtable)

# Average the posterior mean RMSEs for the different seeds
rmse.table <- round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4)
print(rmse.table)

rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xtrue$time) < 1e-6))
for (i in 1:length(ours)) {
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  ours[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
rmse.table <- round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4)
print(rmse.table)

for (i in 1:length(ours)) {
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  ours[[i]]$rmseOdeExpPM <- sqrt(apply((exp(xdesolvePM.obs) - exp(xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
rmse_orig <- round(apply(sapply(ours, function(x) x$rmseOdeExpPM), 1, mean), digits=4)
print(rmse_orig)

rmse_ramsay_log <- round(apply(sapply(ramsayRmseLog, identity), 1, mean, na.rm=TRUE), digits=4)
rmse_ramsay_orig <- round(apply(sapply(ramsayRmseOrig, identity), 1, mean, na.rm=TRUE), digits=4)

xtable(rbind(rmse.table, rmse_ramsay_log))
xtable(rbind(rmse_orig, rmse_ramsay_orig))


logscale <- FALSE
# Make the figures showing Ours using ODE solver results
# use the same axis limits for both methods for easier visual comparison
# on log scale
ylim_lower <- c(-1, -1, -2)
ylim_upper <- c(3, 2, 5)

# on original scale
if(!logscale){
  ylim_lower <- c(0, 0, 0)
  ylim_upper <- c(13, 4, 25)
}

# names of components to use on y-axis label
compnames <- c("P", "M", "H")

desolveOurs <- sapply(ours, function(x) x$xdesolvePM[,-1], simplify = "array")
ourLB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.025))
ourMed <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.5))
ourUB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.975))

xdesolveTRUE <-ours[[1]]$xdesolveTRUE

if(!logscale){
  # on original scale
  # quantile is preserved after exponentiation
  ourLB <- exp(ourLB)
  ourMed <- exp(ourMed)
  ourUB <- exp(ourUB)
  xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
}

times <- ours[[1]]$xdesolveTRUE[,1]
times2 <- seq(0, 240, 0.01)

pdf(width = 20, height = 5, file=paste0(rdaDir, "/plotOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times2, rev(times2)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times2, ourMed[,i])
}
dev.off()


oursExpXdesolvePM <- sapply(oursExpXdesolvePM, function(x) x[,-1], simplify = "array")
ourLB <- apply(oursExpXdesolvePM, c(1,2), function(x) quantile(x, 0.025))
ourMed <- apply(oursExpXdesolvePM, c(1,2), function(x) quantile(x, 0.5))
ourUB <- apply(oursExpXdesolvePM, c(1,2), function(x) quantile(x, 0.975))
xdesolveTRUE <-ours[[1]]$xdesolveTRUE
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
times <- xdesolveTRUE[,1]
pdf(width = 20, height = 5, file=paste0(rdaDir, "/plotOursHes1OrigScale.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times2, rev(times2)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times2, ourMed[,i])
}
dev.off()


# theta posterior mean table 
oursPostTheta <- sapply(oursPostTheta, identity, simplify = "array")

mean_est <- rbind(
  rowMeans(oursPostTheta[,1,])
)

sd_est <- rbind(
  apply(oursPostTheta[,1,], 1, sd)
)

printr <- function(x) format(round(x, 4), nsmall=4)
tablizeEstErr <- function(est, err){
  paste(format(round(est, 4), nsmall=4), "\\pm", format(round(err, 4), nsmall=4))
}

tab <- rbind(
  c("Ours", tablizeEstErr(mean_est[1,],sd_est[1,]))
)
tab <- data.frame(tab)
colnames(tab) <- c("Method", letters[1:7])
rownames(tab) <- NULL
library(xtable)
print(xtable(tab), include.rownames=FALSE)
print(xtable(t(tab)))


# theta posterior credible interval coverage table 
coverage <- rbind(
  printr(rowMeans((oursPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= oursPostTheta[,3,])))
)
coverage <- cbind(c("Ours"), coverage)
colnames(coverage) <- c("Method", letters[1:7])
print(xtable(coverage), include.rownames=FALSE)

print(xtable(cbind(t(tab), t(coverage))))

# X posterior mean plot
oursPostX <- sapply(oursPostX, identity, simplify = "array")

xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
if(!logscale){
  xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
}

pdf(width = 20, height = 5, file=paste0(rdaDir, "/posteriorxOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
  
  if(!logscale){
    # quantile is preserved after exponentiation
    ourLB <- exp(ourLB)
    ourEst <- exp(ourEst)
    ourUB <- exp(ourUB)
  }
  
  times <- xsim$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourEst, col="forestgreen")
}
dev.off()


xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
  
oursPostExpX <- sapply(oursPostExpX, identity, simplify = "array")
pdf(width = 20, height = 5, file=paste0(rdaDir, "/posteriorExpxHes1Ours.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.975)
  
  times <- xsim$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourEst, col="forestgreen")
}
dev.off()

matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)


# smooth visualization with illustration
xdesolveTRUE <-ours[[1]]$xdesolveTRUE
xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
id <- seq(1, nrow(xdesolveTRUE), by=50)
xdesolveTRUE <- xdesolveTRUE[id,]
plot(xdesolveTRUE[,1], xdesolveTRUE[,2], type="l")

ourExpXdesolveLB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.025))
ourExpXdesolveMed <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.5))
ourExpXdesolveUB <- apply(oursExpXdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.975))

phiVisualization <- rbind(
  c(2.07, 0.38, 0.45),
  c(64, 40, 23)
)


ylim_lower <- c(1.5, 0.5, 0)
ylim_upper <- c(9.0, 3.1, 19)

pdf(width = 20, height = 5, file=paste0(rdaDir, "/posteriorExpxHes1OursIllustration.pdf"))
par(mfrow=c(1, ncol(xsim)+1))

matplot(xtrue[, "time"], exp(xtrue[, -1]), type="l", lty=1, col=c(4,6,"goldenrod1"), xlab="time", ylab=NA)
matplot(xsim.obs$time, exp(xsim.obs[,-1]), type="p", col=c(4,6,"goldenrod1"), pch=20, add = TRUE)
mtext('sample observations', cex=1.5)
legend("topright", c("true P", "true M", "true H", "observed P", "observed M"), 
       lty=c(1,1,1,NA,NA), pch=c(NA,NA,NA,20,20), col=c(4,6,"goldenrod1"), cex=1.5)

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


ramsayXdesolvePM <- sapply(ramsayXdesolvePM, identity, simplify = "array")

pdf(width = 20, height = 5, file=paste0(rdaDir, "/posteriorExpxHes1Ramsay.pdf"))
par(mfrow=c(1, ncol(xsim)+1))

matplot(xtrue[, "time"], exp(xtrue[, -1]), type="l", lty=1, col=c(4,6,"goldenrod1"), xlab="time", ylab=NA)
matplot(xsim.obs$time, exp(xsim.obs[,-1]), type="p", col=c(4,6,"goldenrod1"), pch=20, add = TRUE)
mtext('sample observations', cex=1.5)
legend("topright", c("true P", "true M", "true H", "observed P", "observed M"), 
       lty=c(1,1,1,NA,NA), pch=c(NA,NA,NA,20,20), col=c(4,6,"goldenrod1"), cex=1.5)

for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(ramsayXdesolvePM[id,i+1,], 1, quantile, probs = 0.5, na.rm=TRUE)
  ourUB <- apply(ramsayXdesolvePM[id,i+1,], 1, quantile, probs = 0.025, na.rm=TRUE)
  ourLB <- apply(ramsayXdesolvePM[id,i+1,], 1, quantile, probs = 0.975, na.rm=TRUE)
  
  ourEst <- exp(ourEst)
  ourUB <- exp(ourUB)
  ourLB <- exp(ourLB)
  
  times <- xdesolveTRUE[,1]
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i], cex=1.5)
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = "skyblue", lty = 1, density = 10, angle = -45)
  
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourEst, col="forestgreen", lwd=3)
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth",  "median reconstructed trajectory", 
                   "CI on reconstructed trajectory"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "forestgreen", NA), density=c(NA, NA, 40), fill=c(0, 0, "skyblue"),
       border=c(0, 0, "skyblue"), angle=c(NA,NA,-45), x.intersp=c(2.5,2.5, 0),  bty = "n", cex=1.8)
dev.off()
