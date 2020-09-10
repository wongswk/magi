# Summarize the results
library(magi)

config <- list()
config$modelName <- "PTrans"
config$noise <- rep(0.01, 5)

rdaDir <- "../results/ptrans-highnoise/"

seeds <- list.files(rdaDir, pattern='*rda$')
seeds <- unlist(lapply(seeds, function(x) gsub(".*PTrans-([0-9]+).*rda", "\\1", x)))

ptransmodel <- list(
  name= config$modelName,
  fOde=magi:::ptransmodelODE,
  fOdeDx=magi:::ptransmodelDx,
  fOdeDtheta=magi:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(4,6)
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
    rmseOdeMAP <- sqrt(apply((xdesolveMAP.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE)) # compared to true traj
    rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
    
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
ourTime <- c()
for (i in seeds) {
  show(i)
  j <- j+1
  
  load(paste0(rdaDir, config$modelName,"-",i,"-noise", config$noise[1], "-fill0.5.rda"))
  xsim.obs <- xsim[complete.cases(xsim),]
  ours[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  
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
  
  ourTime[j] <- OursTimeUsed 

}

save(ours,oursPostX, oursPostTheta, ourTime, file=paste0(rdaDir,"ourResults.rda"))

# Trajectory RMSE table
rmse.table <- round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4)
show(rmse.table)

# Make the figures for inferred and reconstructed trajectories

# Reconstructed trajectories
showBlackMedian <- FALSE
ylim_lower <- c(0,0,0.3,0,0)
ylim_upper <- c(1, 0.5, 1.1, 0.4, 0.65)
# names of components to use on y-axis label
compnames <- c("S", "Sd", "R", "RS", "Rpp")

desolveOurs <- sapply(ours, function(x) x$xdesolvePM[,-1], simplify = "array")
ourLB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.025))
ourMed <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.5))
ourUB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.975))

times <- ours[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(rdaDir, "plotOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, ours[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]), cex.axis=1.2)
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, ours[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  if(showBlackMedian){
    lines(times, ourMed[,i])
  }
  mtext(compnames[i], cex=2)
}
dev.off()


# Inferred trajectories
oursPostX <- simplify2array(oursPostX)
xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(rdaDir, "posteriorxOurs.pdf"))
layout(rbind(c(1,2,3,4,5), c(6,6,6,6,6)), heights = c(8,1))

  i <- 1 
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

  zoomin_x <- 2*c(37, 43)
  zoomin_y <- c(-0.012, 0.012)
  polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
          col = NA, border = 1)
  
  zoomout_x <- 2*c(25, 80)
  zoomout_y <- c(0.25, 0.55)
  polygon(c(times[zoomout_x], rev(times[zoomout_x])), rep(zoomout_y, each=2),
          col = NA, border = 1)
  
  lines(times[c(zoomin_x[1], zoomout_x[1])], c(zoomin_y[2], zoomout_y[1]))
  lines(times[c(zoomin_x[2], zoomout_x[2])], c(zoomin_y[2], zoomout_y[1]))
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
  
    
  i <- 2 
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
  
  i <- 3 
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

  zoomin_x <- 2*c(55, 65)
  zoomin_y <- c(0.73, 0.83)
  polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
          col = NA, border = 1)
  
  zoomout_x <- 2*c(50, 86)
  zoomout_y <- c(0.3, 0.6)
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
    
  i <- 4 
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
  
  zoomin_x <- 2*c(35, 45)
  zoomin_y <- c(-0.012, 0.012)
  polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
          col = NA, border = 1)
  
  zoomout_x <- 2*c(25, 80)
  zoomout_y <- c(0.1, 0.22)
  polygon(c(times[zoomout_x], rev(times[zoomout_x])), rep(zoomout_y, each=2),
          col = NA, border = 1)
  
  lines(times[c(zoomin_x[1], zoomout_x[1])], c(zoomin_y[2], zoomout_y[1]))
  lines(times[c(zoomin_x[2], zoomout_x[2])], c(zoomin_y[2], zoomout_y[1]))
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
  
  i <- 5
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
  
  zoomin_x <- 2*c(50, 60)
  zoomin_y <- c(0.21, 0.312)
  polygon(c(times[zoomin_x], rev(times[zoomin_x])), rep(zoomin_y, each=2),
          col = NA, border = 1)
  
  zoomout_x <- 2*c(60, 90)
  zoomout_y <- c(0.37, 0.65)
  polygon(c(times[zoomout_x], rev(times[zoomout_x])), rep(zoomout_y, each=2),
          col = NA, border = 1)
  
  lines(times[c(zoomin_x[1], zoomout_x[1])], c(zoomin_y[2], zoomout_y[1]))
  lines(times[c(zoomin_x[2], zoomout_x[2])], c(zoomin_y[2], zoomout_y[1]))
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
  xcoords <- c(0, 0.5, 1)
  secondvector <- (1:length(xcoords))-1
  textwidths <- xcoords/secondvector # this works for all but the first element
  textwidths[1] <- 0
  
  legend("center", c("truth", "median of all inferred trajectories", "95% interval from the 2.5 and 97.5 percentile of all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
         col = c("red", "forestgreen", NA), fill=c(0, 0,"skyblue"), text.width=c(0, 0.4, 0.05), bty = "n",
         border=c(0, 0, "skyblue"), pch=c(NA, NA, 15), cex=1.8, horiz=TRUE)
  dev.off()


