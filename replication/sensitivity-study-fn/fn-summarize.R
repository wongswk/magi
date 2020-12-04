# Summarize the results
library(magi)

config <- list()
config$modelName <- "FN"
config$noise <- rep(0.2, 2)
rdaDir <- "../results/fn-fill2-nobs41/"

sink(paste0(rdaDir, "/result.txt"))

seeds <- list.files(rdaDir)  ## get the list of seeds ran
seeds <- seeds[grep(".*FN-([0-9]+)-noise.*", seeds)]
seeds <- unique(gsub(".*FN-([0-9]+)-noise.*", "\\1", seeds))

fnmodel <- list(
  fOde=magi:::fODE,
  fOdeDx=magi:::fnmodelDx,
  fOdeDtheta=magi:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

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
j <- 0
ours <- list()

oursPostX <- list()

oursPostTheta <- list()

for (i in seeds) {
  show(i)
  j <- j+1
  
  load(paste0(rdaDir, config$modelName,"-",i,"-noise", config$noise[1], ".rda"))
  xsim.obs <- xsim[complete.cases(xsim),]
  ours[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
  oursPostX[[j]] <- cbind(
    apply(gpode$xsampled, 2:3, mean), 
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  oursPostTheta[[j]] <- cbind(
    apply(gpode$theta, 2, mean), 
    apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  )
  
}

####### Recalculate RMSE against true trajectory (load compare.rda first)
rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xtrue$time) < 1e-6))
for (i in 1:length(ours)) {
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  ours[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
}

save.image(paste0(rdaDir, "fn-ours-all.rda"))

# Average the posterior mean RMSEs for the different seeds
rmse.table <- rbind( round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4))
print(rmse.table)

# Make the figures comparing Wenk and Ours using ODE solver results
# use the same axis limits for both methods for easier visual comparison
ylim_lower <- c(-2.5,-1.5)
ylim_upper <- c(2.5,1.5)
# names of components to use on y-axis label
compnames <- c("V", "R")

desolveOurs <- sapply(ours, function(x) x$xdesolvePM[,-1], simplify = "array")
ourLB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.025))
ourMed <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.5))
ourUB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.975))

times <- ours[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(rdaDir, "plotOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, ours[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, ours[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourMed[,i])
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

rmse_est <- rbind(
  apply(oursPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2)))
)


digit_precision = 2
printr <- function(x, precision=digit_precision) format(round(x, precision), nsmall=precision)
tablizeEstErr <- function(est, err){
  paste(format(round(est, digit_precision), nsmall=digit_precision), "\\pm", format(round(err, digit_precision), nsmall=digit_precision))
}

tab <- rbind(
  c("Ours", tablizeEstErr(mean_est[1,],sd_est[1,]))
)
tab <- data.frame(tab)
colnames(tab) <- c("Method", "a", "b", "c")
rownames(tab) <- NULL
library(xtable)
print(xtable(tab), include.rownames=FALSE)



# theta posterior credible interval coverage table 
coverage <- rbind(
  printr(rowMeans((oursPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= oursPostTheta[,3,])))
)
coverage <- cbind(c("Ours"), coverage)
colnames(coverage) <- c("Method", "a", "b", "c")
print(xtable(coverage), include.rownames=FALSE)

rmse_theta <- rbind(
  printr(apply(oursPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))), 4)
)
rmse_theta <- cbind(c("Ours"), rmse_theta)
colnames(rmse_theta) <- c("Method", "a", "b", "c")
print(xtable(rmse_theta), include.rownames=FALSE)

# X posterior mean plot
oursPostX <- sapply(oursPostX, identity, simplify = "array")

xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)

pdf(width = 20, height = 5, file=paste0(rdaDir, "posteriorxOurs.pdf"))
layout(t(1:3), widths = c(2,2,1))
i = 1
ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
times <- xsim$time

plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
mtext(compnames[i], cex=2)
polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
        col = "skyblue", border = NA)
lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
lines(times, ourEst, col="forestgreen", lwd=3)

zoomin_x <- c(85, 89)
zoomin_y <- c(1.80, 2.05)
polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
        col = NA, border = 1)

zoomout_x <- c(80, 116)
zoomout_y <- c(-2.5, -0.4)
polygon(c(times[zoomout_x], rev(times[zoomout_x])), rep(zoomout_y, each=2),
        col = NA, border = 1)

lines(times[c(zoomin_x[1], zoomout_x[1])], c(zoomin_y[1], zoomout_y[2]))
lines(times[c(zoomin_x[2], zoomout_x[2])], c(zoomin_y[1], zoomout_y[2]))
zoomout_id <- zoomout_x[1]:zoomout_x[2]
zoomin_id <- zoomin_x[1]:zoomin_x[2]
zoomtrans_y <- function(y) {
  y <- y[zoomin_id]
  (y - zoomin_y[1])/diff(zoomin_y) * diff(zoomout_y) + zoomout_y[1]
}
zoomtrans_x <- function(x) {
  zoomin_x <- times[zoomin_x]
  zoomout_x <- times[zoomout_x]
  (x - zoomin_x[1])/diff(zoomin_x) * diff(zoomout_x) + zoomout_x[1]
}
polygon(c(zoomtrans_x(times[zoomin_id]), rev(zoomtrans_x(times[zoomin_id]))), c(zoomtrans_y(ourUB), rev(zoomtrans_y(ourLB))),
        col = "skyblue", border = NA)
lines(zoomtrans_x(times[zoomin_id]), zoomtrans_y(xdesolveTRUE[,1+i]), col="red", lwd=4)
lines(zoomtrans_x(times[zoomin_id]), zoomtrans_y(ourEst), col="forestgreen", lwd=3)

i = 2
ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
times <- xsim$time

plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
mtext(compnames[i], cex=2)
polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
        col = "skyblue", border = NA)
lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
lines(times, ourEst, col="forestgreen", lwd=3)

zoomin_x <- c(73, 79)
zoomin_y <- c(0.9, 1.10)
polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
        col = NA, border = 1)

zoomout_x <- c(57, 94)
zoomout_y <- c(-1.5, -0.4)
polygon(c(times[zoomout_x], rev(times[zoomout_x])), rep(zoomout_y, each=2),
        col = NA, border = 1)

lines(times[c(zoomin_x[1], zoomout_x[1])], c(zoomin_y[1], zoomout_y[2]))
lines(times[c(zoomin_x[2], zoomout_x[2])], c(zoomin_y[1], zoomout_y[2]))
zoomout_id <- zoomout_x[1]:zoomout_x[2]
zoomin_id <- zoomin_x[1]:zoomin_x[2]
zoomtrans_y <- function(y) {
  y <- y[zoomin_id]
  (y - zoomin_y[1])/diff(zoomin_y) * diff(zoomout_y) + zoomout_y[1]
}
zoomtrans_x <- function(x) {
  zoomin_x <- times[zoomin_x]
  zoomout_x <- times[zoomout_x]
  (x - zoomin_x[1])/diff(zoomin_x) * diff(zoomout_x) + zoomout_x[1]
}
polygon(c(zoomtrans_x(times[zoomin_id]), rev(zoomtrans_x(times[zoomin_id]))), c(zoomtrans_y(ourUB), rev(zoomtrans_y(ourLB))),
        col = "skyblue", border = NA)
lines(zoomtrans_x(times[zoomin_id]), zoomtrans_y(xdesolveTRUE[,1+i]), col="red", lwd=4)
lines(zoomtrans_x(times[zoomin_id]), zoomtrans_y(ourEst), col="forestgreen", lwd=3)


par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth", "median of all inferred trajectories", "2.5% to 97.5% percentile\nof all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "forestgreen", NA), fill=c(0, 0, "skyblue"), pch=c(NA, NA, 15), x.intersp=c(2.5,2.5,0),
       border=c(0, 0, "skyblue"), bty = "n", cex=1.8)
dev.off()

sink()
