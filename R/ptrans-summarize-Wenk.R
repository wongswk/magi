# Summarize the results
library(gpds)

config <- list()
config$modelName <- "PTrans"
config$noise <- rep(0.01, 5)
rdaDir <- "../DynSysResults/PTrans-noise0.001/results/"   ## where ours & Dondel rda saved
outDirWenk <- "../DynSysResults/PTrans-noise0.001/"   ## where all the seeds are in separate folders with Wenk's output

#rdaDir <- "../results/cpp/"
#outDirWenk <- rdaDir

seeds <- list.files(outDirWenk, pattern='^\\d+$')  ## get the list of seeds ran

#seeds <- list.files(outDirWenk, pattern='*rda$')
seeds <- unlist(lapply(seeds, function(x) gsub(".*PTrans-([0-9]+).*rda", "\\1", x)))

ptransmodel <- list(
  name= config$modelName,
  fOde=gpds:::ptransmodelODE,
  fOdeDx=gpds:::ptransmodelDx,
  fOdeDtheta=gpds:::ptransmodelDtheta,
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
Wenk <- list()
WenkPostX <- list()
Dondel <- list()
for (i in seeds) {
  show(i)
  j <- j+1
  
  load(paste0(rdaDir, config$modelName,"-",i,"-noise", config$noise[1], ".rda"))
  xsim.obs <- xsim[complete.cases(xsim),]
  # ours[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  # 
  # oursPostX[[j]] <- cbind(
  #   apply(gpode$xsampled, 2:3, mean), 
  #   apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
  #   apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  # )
  # oursPostTheta[[j]] <- cbind(
  #   apply(gpode$theta, 2, mean), 
  #   apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
  #   apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  # )
  
  
  # load(paste0(rdaDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".rda"))
  # config$n.iter.Dondel <- config$n.iter.Dondel / 25
  # x_means <- apply(xsim.obs[,-1],2,mean)
  # 
  # xsampled <- (sapply(agm.result$x.samples, function(x) x, simplify="array"))
  # xsampled <- aperm(apply(xsampled, c(1,2), function(x) x + x_means), c(2,3,1))
  # thetasampled <- agm.result$posterior.samples
  # lglik <- agm.result$ll
  # 
  # burnin <- as.integer(config$n.iter.Dondel*config$burninRatio)
  # 
  # gpode <- list(theta= thetasampled[-(1:burnin),],
  #               xsampled= xsampled[-(1:burnin),,],
  #               lglik=  lglik[-(1:burnin)],
  #               sigma=  matrix(0, nrow = config$n.iter-burnin, ncol = ncol(xsim)-1)) # not sampled in this method
  # gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  #   with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  # gpode$fode <- aperm(gpode$fode, c(3,1,2))
  # Dondel[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  # 
  dataDir <- paste0(outDirWenk, i, "/")

  samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
  xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
  xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
  thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))
  lliklist <- as.vector(as.matrix(read.table(paste0(dataDir, "lliklist.csv") )))

  llikId <- 0  ### llik is in its own file
  xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim.obs[,-1])))
  thetaId <- (max(xId)+1):(max(xId)+length(ptransmodel$thetaLowerBound))
  #sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))  ## no sigma sampled

  ## Un-standardize their X's and untransform thetas
  samplesCpp[,xId] <- t(apply(samplesCpp[,xId],1, function(x) xSD*x + xMean))
  samplesCpp[,thetaId] <- t(apply(samplesCpp[,thetaId],1, function(x) x * 10^thetaMag))

  burnin <- as.integer(config$n.iter.Wenk*config$burninRatio)
  gpode <- list(theta= samplesCpp[-(1:burnin), thetaId],
                xsampled=array(samplesCpp[-(1:burnin), xId],
                               dim=c(nrow(samplesCpp)-burnin, nrow(xsim.obs), ncol(xsim)-1)),
                lglik=  lliklist[-(1:burnin)], ###samplesCpp[llikId,-(1:burnin)],
                sigma=  matrix(0, nrow = nrow(samplesCpp) - burnin, ncol = ncol(xsim)-1)) # t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
  gpode$fode <- sapply(1:length(gpode$lglik), function(t)
    with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  gpode$fode <- aperm(gpode$fode, c(3,1,2))

  Wenk[[j]]<- rmsePostSamples( xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  
  WenkPostX[[j]] <- cbind(
    apply(gpode$xsampled, 2:3, mean), 
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  

}

save(Wenk, WenkPostX, file=paste0(outDirWenk,"Wenk.rda"))
#save(ours,oursPostX, oursPostTheta, Wenk,Dondel, file=paste0(outDirWenk,"compare.rda"))

####### Recalculate RMSE against true trajectory (load compare.rda first)
# rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xtrue$time) < 1e-6))
# for (i in 1:length(ours)) {
#     xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
#     xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
#     ours[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
# }
# for (i in 1:length(Wenk)) {
#   xdesolveTRUE.obs <- Wenk[[i]]$xdesolveTRUE[rowId,-1]
#   xdesolvePM.obs <- Wenk[[i]]$xdesolvePM[rowId,-1]
#   Wenk[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
# }
# for (i in 1:length(Dondel)) {
#   xdesolveTRUE.obs <- Dondel[[i]]$xdesolveTRUE[rowId,-1]
#   xdesolvePM.obs <- Dondel[[i]]$xdesolvePM[rowId,-1]
#   Dondel[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
# }
##########  end recalculation

# Average the posterior mean RMSEs for the different seeds
#rmse.table <- rbind( round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4),
#  round(apply(sapply(Wenk, function(x) x$rmseOdePM), 1, mean), digits=4),
#  round(apply(sapply(Dondel, function(x) x$rmseOdePM), 1, mean), digits=4))

rmse.table <- round(apply(sapply(Wenk, function(x) x$rmseOdePM), 1, mean), digits=4)

# Make the figures comparing Wenk and Ours using ODE solver results
# use the same axis limits for both methods for easier visual comparison
showBlackMedian <- FALSE

ylim_lower <- c(0,0,0.3,0,0)
ylim_upper <- c(1, 0.5, 1.1, 0.4, 0.65)
# names of components to use on y-axis label
compnames <- c("S", "Sd", "R", "RS", "Rpp")

# desolveOurs <- sapply(ours, function(x) x$xdesolvePM[,-1], simplify = "array")
# ourLB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.025))
# ourMed <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.5))
# ourUB <- apply(desolveOurs, c(1,2), function(x) quantile(x, 0.975))
# 
# times <- ours[[1]]$xdesolveTRUE[,1]
# 
# pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotOurs.pdf"))
# par(mfrow=c(1, ncol(xsim)-1))
# for (i in 1:(ncol(xsim)-1)) {
#   plot(times, ours[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
#   polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
#           col = "grey80", border = NA)
#   lines(times, ours[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
#   lines(times, ourMed[,i])
# }
# dev.off()

desolveWenk <- sapply(Wenk, function(x) x$xdesolvePM[,-1], simplify = "array")
WenkLB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.025))
WenkMed <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.5))
WenkUB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.975))

times <- Wenk[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotWenk.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, Wenk[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]), cex.axis=1.2)
  polygon(c(times, rev(times)), c(WenkUB[,i], rev(WenkLB[,i])),
          col = "grey80", border = NA)
  lines(times, Wenk[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  if(showBlackMedian){
    lines(times, WenkMed[,i])
  }
  mtext(compnames[i], cex=2)
}
dev.off()

# X posterior mean plot
#oursPostX <- sapply(oursPostX, identity, simplify = "array")
oursPostX <- simplify2array(WenkPostX)

# use the same axis limits for both methods for easier visual comparison
ylim_lower <- c(0,0,0.3,0,0)
ylim_upper <- c(1, 0.25, 1, 0.4, 0.65)
# names of components to use on y-axis label
compnames <- c("S", "Sd", "R", "RS", "Rpp")


xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxWenk.pdf"))
#layout(t(1:6), widths = c(2,2,2,2,2,0.5))
layout(t(1:5), widths = c(2,2,2,2,2))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
  times <- xsim.obs$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[i], cex=2)
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(xsim$time, xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourEst, col="forestgreen", lwd=3)
}
# par(mar=rep(0,4))
# plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
# legend("center", c("truth", "median posterior mean", "95% interval"), lty=c(1,1,0), lwd=c(4,3,0),
#        col = c("red", "forestgreen", NA), fill=c(0, 0, "skyblue"), pch=c(NA, NA, 15),
#        border=c(0, 0, "skyblue"), bty = "n", cex=1.8)
dev.off()


