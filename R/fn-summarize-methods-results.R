# Summarize the results
library(gpds)

config <- list()
config$modelName <- "FN"
config$noise <- rep(0.2, 2)
rdaDir <- "comparison/results/"   ## where ours & Dondel rda saved
outDirWenk <- "comparison/results/"   ## where all the seeds are in separate folders with Wenk's output

seeds <- list.files(outDirWenk, pattern='^\\d+$')  ## get the list of seeds ran

fnmodel <- list(
  fOde=gpds:::fODE,
  fOdeDx=gpds:::fnmodelDx,
  fOdeDtheta=gpds:::fnmodelDtheta,
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
Wenk <- list()
Dondel <- list()
Ramsay <- list()

oursPostX <- list()
WenkPostX <- list()
DondelPostX <- list()
RamsayPostX <- list()

oursPostTheta <- list()
WenkPostTheta <- list()
DondelPostTheta <- list()
RamsayPostTheta <- list()

for (i in seeds) {
  if(!file.exists(paste0(outDirWenk, i, "/MCMCMatrix.csv"))){
    next
  }
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
  
  load(paste0(rdaDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".rda"))
  config$n.iter.Dondel <- config$n.iter.Dondel / 25
  x_means <- apply(xsim.obs[,-1],2,mean)
  
  xsampled <- (sapply(agm.result$x.samples, function(x) x, simplify="array"))
  xsampled <- aperm(apply(xsampled, c(1,2), function(x) x + x_means), c(2,3,1))
  thetasampled <- agm.result$posterior.samples
  lglik <- agm.result$ll
  
  burnin <- as.integer(config$n.iter.Dondel*config$burninRatio)
  
  gpode <- list(theta= thetasampled[-(1:burnin),],
                xsampled= xsampled[-(1:burnin),,],
                lglik=  lglik[-(1:burnin)],
                sigma=  matrix(0, nrow = config$n.iter-burnin, ncol = ncol(xsim)-1)) # not sampled in this method
  gpode$fode <- sapply(1:length(gpode$lglik), function(t)
    with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  gpode$fode <- aperm(gpode$fode, c(3,1,2))
  Dondel[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  DondelPostX[[j]] <- cbind(
    apply(gpode$xsampled, 2:3, mean),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  DondelPostTheta[[j]] <- cbind(
    apply(gpode$theta, 2, mean),
    apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  )
  
  dataDir <- paste0(outDirWenk, config$seed, "/")
  
  samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
  xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
  xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
  thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))
  lliklist <- as.vector(as.matrix(read.table(paste0(dataDir, "lliklist.csv") )))
  
  llikId <- 0  ### llik is in its own file
  xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim.obs[,-1])))
  thetaId <- (max(xId)+1):(max(xId)+length(fnmodel$thetaLowerBound))
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
    with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  gpode$fode <- aperm(gpode$fode, c(3,1,2))
  
  
  Wenk[[j]]<- rmsePostSamples( xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
  WenkPostX[[j]] <- cbind(
    apply(gpode$xsampled, 2:3, mean),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  WenkPostTheta[[j]] <- cbind(
    apply(gpode$theta, 2, mean),
    apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  )
  
  gpodeRamsay <- readRDS(paste0(rdaDir, config$modelName,"-",i,"-noise", config$noise[1], "-gpodeRamsay.rds"))
  gpodeRamsay$sigma <- matrix(pram.true$sigma, nrow=length(gpodeRamsay$lglik), ncol=length(pram.true$sigma), byrow = TRUE)
  xsimRamsay <- data.frame(time=seq(0, 20, 0.1), V=NA, R=NA)
  xsimRamsay$V[match(na.omit(xsim)$time, xsimRamsay$time)] <- na.omit(xsim)[,2]
  xsimRamsay$R[match(na.omit(xsim)$time, xsimRamsay$time)] <- na.omit(xsim)[,3]
  gpodeRamsay$xsampled <- gpodeRamsay$xsampled[,2:202,]
  gpodeRamsay$fode <- sapply(1:length(gpodeRamsay$lglik), function(t)
    with(gpodeRamsay, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
  gpodeRamsay$fode <- aperm(gpodeRamsay$fode, c(3,1,2))
  
  Ramsay[[j]] <- rmsePostSamples(xtrue, dotxtrue, xsim, gpodeRamsay, pram.true, config, odemodel)
  RamsayPostX[[j]] <- cbind(
    apply(gpodeRamsay$xsampled, 2:3, mean), 
    apply(gpodeRamsay$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpodeRamsay$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  RamsayPostTheta[[j]] <- cbind(
    apply(gpodeRamsay$theta, 2, mean), 
    apply(gpodeRamsay$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpodeRamsay$theta, 2, function(x) quantile(x, 0.975))
  )
}

save.image(file=paste0(outDirWenk,"compare-withRamsay.rda"))

####### Recalculate RMSE against true trajectory (load compare.rda first)
rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xtrue$time) < 1e-6))
for (i in 1:length(ours)) {
  xdesolveTRUE.obs <- ours[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
  ours[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
for (i in 1:length(Wenk)) {
  xdesolveTRUE.obs <- Wenk[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- Wenk[[i]]$xdesolvePM[rowId,-1]
  Wenk[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
}
for (i in 1:length(Dondel)) {
  xdesolveTRUE.obs <- Dondel[[i]]$xdesolveTRUE[rowId,-1]
  xdesolvePM.obs <- Dondel[[i]]$xdesolvePM[rowId,-1]
  Dondel[[i]]$rmseOdePM <- sqrt(apply((xdesolvePM.obs - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))   # compared to true traj
}


# Average the posterior mean RMSEs for the different seeds
rmse.table <- rbind( round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=3),
                     round(apply(sapply(Ramsay, function(x) x$rmseOdePM), 1, mean), digits=3),
                     round(apply(sapply(Wenk, function(x) x$rmseOdePM), 1, mean), digits=3),
                     round(apply(sapply(Dondel, function(x) x$rmseOdePM), 1, mean), digits=3))
print("trajectory RMSE")
print(cbind(c("ours", "ramsay", "wenk", "dondel"), rmse.table))

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

pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, ours[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, ours[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourMed[,i])
}
dev.off()


desolveRamsay <- sapply(Ramsay, function(x) x$xdesolvePM[,-1], simplify = "array")
RamsayLB <- apply(desolveRamsay, c(1,2), function(x) quantile(x, 0.025))
RamsayMed <- apply(desolveRamsay, c(1,2), function(x) quantile(x, 0.5))
RamsayUB <- apply(desolveRamsay, c(1,2), function(x) quantile(x, 0.975))

times <- Ramsay[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotRamsay.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, Ramsay[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(RamsayUB[,i], rev(RamsayLB[,i])),
          col = "grey80", border = NA)
  lines(times, Ramsay[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, RamsayMed[,i])
}
dev.off()


desolveWenk <- sapply(Wenk, function(x) x$xdesolvePM[,-1], simplify = "array")
WenkLB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.025))
WenkMed <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.5))
WenkUB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.975))

times <- Wenk[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotWenk.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, Wenk[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(WenkUB[,i], rev(WenkLB[,i])),
          col = "grey80", border = NA)
  lines(times, Wenk[[1]]$xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, WenkMed[,i])
}
dev.off()


# theta posterior mean table 
oursPostTheta <- sapply(oursPostTheta, identity, simplify = "array")
RamsayPostTheta <- sapply(RamsayPostTheta, identity, simplify = "array")
WenkPostTheta <- sapply(WenkPostTheta, identity, simplify = "array")
DondelPostTheta <- sapply(DondelPostTheta, identity, simplify = "array")

mean_est <- rbind(
  rowMeans(oursPostTheta[,1,]),
  rowMeans(WenkPostTheta[,1,]),
  rowMeans(DondelPostTheta[,1,]),
  rowMeans(RamsayPostTheta[,1,])
)

sd_est <- rbind(
  apply(oursPostTheta[,1,], 1, sd),
  apply(WenkPostTheta[,1,], 1, sd),
  apply(DondelPostTheta[,1,], 1, sd),
  apply(RamsayPostTheta[,1,], 1, sd)
)

rmse_est <- rbind(
  apply(oursPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))),
  apply(WenkPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))),
  apply(DondelPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))),
  apply(RamsayPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2)))
)


digit_precision = 2
printr <- function(x, precision=digit_precision) format(round(x, precision), nsmall=precision)
tablizeEstErr <- function(est, err){
  paste(format(round(est, digit_precision), nsmall=digit_precision), "\\pm", format(round(err, digit_precision), nsmall=digit_precision))
}

tab <- rbind(
  c("Ours", tablizeEstErr(mean_est[1,],sd_est[1,])),
  c("Wenk", tablizeEstErr(mean_est[2,],sd_est[2,])),
  c("Dondel", tablizeEstErr(mean_est[3,],sd_est[3,]))
)
tab <- data.frame(tab)
colnames(tab) <- c("Method", "a", "b", "c")
rownames(tab) <- NULL
library(xtable)
print(xtable(tab), include.rownames=FALSE)



# theta posterior credible interval coverage table 
coverage <- rbind(
  printr(rowMeans((oursPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= oursPostTheta[,3,]))),
  printr(rowMeans((WenkPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= WenkPostTheta[,3,]))),
  printr(rowMeans((DondelPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= DondelPostTheta[,3,])))
)
coverage <- cbind(c("Ours", "Wenk", "Dondel"), coverage)
colnames(coverage) <- c("Method", "a", "b", "c")
print(xtable(coverage), include.rownames=FALSE)

rmse_theta <- rbind(
  printr(apply(oursPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))), 3),
  printr(apply(WenkPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))), 3),
  printr(apply(DondelPostTheta[,1,] - pram.true$theta, 1, function(x) sqrt(mean(x^2))), 3)
)
rmse_theta <- cbind(c("Ours", "Wenk", "Dondel"), rmse_theta)
colnames(rmse_theta) <- c("Method", "a", "b", "c")
print(xtable(rmse_theta), include.rownames=FALSE)


inference <- cbind(tab[,1,drop=FALSE],cbind(tab[,-1], coverage[,-1])[,c(1,4,2,5,3,6)])
print(xtable(inference), include.rownames=FALSE)

inference <- cbind(t(tab[,-1]), t(coverage[,-1]))[,c(1,4,2,5,3,6)]
print(xtable(inference))

print("main table below")
inference <- cbind(t(tab[,-1]), t(rmse_theta[,-1]))[,c(1,4,2,5,3,6)]
print(xtable(inference))

inference <- cbind(t(tab[,-1]), t(rmse_theta[,-1]), t(coverage[,-1]))[,c(1,4,7,2,5,8,3,6,9)]
print(xtable(inference))

# X posterior mean plot
oursPostX <- sapply(oursPostX, identity, simplify = "array")
WenkPostX <- sapply(WenkPostX, identity, simplify = "array")
DondelPostX <- sapply(DondelPostX, identity, simplify = "array")

xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxOurs.pdf"))
layout(t(1:3), widths = c(2,2,1))
for (i in 1:(ncol(xsim)-1)) {
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
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth", "median posterior mean", "95% interval on posterior mean"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "forestgreen", NA), fill=c(0, 0, "skyblue"), pch=c(NA, NA, 15),
       border=c(0, 0, "skyblue"), bty = "n", cex=1.8)
dev.off()

xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim.obs$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxWenk.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(WenkPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(WenkPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(WenkPostX[,i,], 1, quantile, probs = 0.975)
  times <- xsim.obs$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourEst, col="forestgreen")
}
dev.off()

