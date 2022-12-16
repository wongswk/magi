# Summarize the results
library(magi)

expt <- 'noisemin0.05-nobs35-uneven'
rdaDir <- paste0("../results/lac-operon/", expt, "/")
pdf_files <- list.files(rdaDir)
rda_files <- pdf_files[grep("lac-operon.*\\.rda", pdf_files)]
pdf_files <- pdf_files[grep("lac-operon.*\\.pdf", pdf_files)]
  
print(length(rda_files))
  
config <- list()
config$modelName <- "lac-operon"
#config$noise <- rep(0.2, 6)

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
oursRecon <- list()
ourTime <- c()

jj <- 0
for (i in rda_files) {
  show(i)
  jj <- jj+1

  load(paste0(rdaDir, i))
  xsim.obs <- xsim[complete.cases(xsim[,1:3]),]
  ours[[jj]] <- rmsePostSamples(xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)

  oursPostX[[jj]] <- cbind(
    apply(gpode$xsampled, 2:3, mean),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
    apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975))
  )
  oursPostTheta[[jj]] <- cbind(
    apply(gpode$theta, 2, mean),
    apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
    apply(gpode$theta, 2, function(x) quantile(x, 0.975))
  )

  oursRecon[[jj]] <- recon[,-1]

  ourTime[jj] <- OursTimeUsed


}

print(paste0(rdaDir,"summary-lacoperon.rda"))
save.image(file=paste0(rdaDir,"summary-lacoperon.rda"))

#load(file=paste0(rdaDir,"summary-lacoperon.rda"))

jj <- 1 #sample(1:100, 1)
load(paste0(rdaDir, rda_files[jj]))
 
xtrue <- xtrue[xtrue$time>=1,]
xtrue <- xtrue[seq(1,nrow(xtrue), by=100),]

pdf(width = 12, height = 6, file=paste0(rdaDir, "lac-inferred-example.pdf"))
compnames <- c("r_I", "I", "Lactose", "ILactose", "Op", "IOp", "RNAP", "RNAPo", "r", "Z")
par(oma=c(0,1.5,0,0))
yaxisLims <- matrix(NA, nrow=2, ncol = 10)
layout(rbind(c(1:5), c(6:10), c(11,11,11,11,11)), heights = c(8,8,1))
for (ii in 1:10) {
  
  par(mar = c(4, 2.5, 1.75, 0.1))
  ourEst <- oursPostX[[jj]][,ii]
  ourEst <- magi:::getMeanCurve(xsim$time, ourEst, xtrue[,1],
                                    t(phiExogenous[,ii]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE)
  
  ourUB <- oursPostX[[jj]][,20+ii]
  ourUB <- magi:::getMeanCurve(xsim$time, ourUB, xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE)
  
  ourLB <- oursPostX[[jj]][,10+ii]
  ourLB <- magi:::getMeanCurve(xsim$time, ourLB, xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE)
  #if (ii == 7)
  #  plot( c(min(xtrue$time),max(xtrue$time)), c(min(c(ourLB,xsim.obs[,ii+1])), max(c(xsim.obs[,ii+1],ourUB))), type='n', xlab='', ylab='', ylim=c(88,108))
  #else  
  yaxisLims[,ii] <- c(min(c(ourLB,xsim.obs[,ii+1])), max(c(xsim.obs[,ii+1],ourUB)))
  if (ii == 7)
    yaxisLims[,ii] <- c(68,128)
  
  plot( c(min(xtrue$time),max(xtrue$time)), yaxisLims[,ii], type='n', xlab='', ylab='')
  
  polygon(c(xtrue[,1], rev(xtrue[,1])), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)    
  # polygon(c(xsim[,1], rev(xsim[,1])), c(ourUB, rev(ourLB)),
  #        col = "skyblue", border = NA)    
  
  #if (ii %in% c(1,6))
  #  title(ylab='Level', cex.lab = 1.5)
  
  if (ii == 8)
    title(xlab='Time (sec)', cex.lab = 1.5)
  
  lines(xtrue[, "time"], xtrue[,ii+1],col='red', lwd=2)
  lines(xtrue[,1], ourEst, col='forestgreen', lwd=1.5)
  mtext(compnames[ii], cex=1.25)
  points(xsim.obs$time, xsim.obs[,ii+1], col='black', pch=16)
  
  if (ii ==1) mtext("       Concentration (arb. unit)",side=2,line=0,outer=TRUE)
}

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n',
     xlab = NA, ylab = NA, frame.plot = FALSE)

legend("center", c("truth", "inferred trajectory",
                   "95% interval", "noisy observations"),
       lty = c(1, 1, 0, 0), lwd = c(2, 2, 0, 1), bty = "n",
       col = c("red", "forestgreen", NA, "black"), fill = c(0, 0, "skyblue", 0),
       border = c(0, 0, "skyblue", 0), pch = c(NA, NA, 15, 16), horiz = TRUE, cex=1.25)


dev.off()



# Inferred trajectories
# oursPostExpX <- simplify2array(oursPostExpX)
# xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = c(0,xsim$time), func = odemodel$modelODE, parms = pram.true$theta)
# xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
# pdf(width = 17, height = 10, file=paste0(rdaDir, "posteriorxOurs.pdf"))
# 
# layout(rbind(c(1,2,3), c(4,5,6), c(7,7,7)), heights = c(8,8,1))
# compnames <- c("m_laci", "m_tetr", "m_ci", "p_laci (unobserved)", "p_tetr (unobserved)", "p_ci (unobserved)")
# for (i in 1:6) {
#   ourEst <- apply(oursPostExpX[,i,], 1, function(x) quantile(x, 0.5))
#   ourUB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.025)
#   ourLB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.975)
#   times <- xsim$time
#   
#   plot(times, ourEst, type="n", xlab="time", ylab=compnames[i]) #, ylim=c(ylim_lower[i], ylim_upper[i]))
#   mtext(compnames[i], cex=1.25)
#   polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
#           col = "skyblue", border = NA)
#   lines(times, xdesolveTRUE[-1,1+i], col="red", lwd=2)
#   lines(times, ourEst, col="forestgreen", lwd=1)
# }
# 
# par(mar=rep(0,4))
# plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
# legend("center", c("truth", "median of all inferred trajectories", "95% interval from the 2.5 and 97.5 percentile of all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
#        col = c("red", "forestgreen", NA), fill=c(0, 0,"skyblue"), text.width=c(0, 0.15, 0.01), bty = "n",
#        border=c(0, 0, "skyblue"), pch=c(NA, NA, 15), cex=1.2, horiz=TRUE)
# 
# dev.off()

#load(file=paste0(rdaDir,"lac-operon_summary.rda"))

# Reconstructed trajectories
#compnames <- paste("Component", 1:10)
oursRecon <- sapply(oursRecon, identity, simplify = "array")

ourLB <- apply(oursRecon, c(1,2), function(x) quantile(x, 0.025))
ourMed <- apply(oursRecon, c(1,2), function(x) quantile(x, 0.5))
ourUB <- apply(oursRecon, c(1,2), function(x) quantile(x, 0.975))

pdf(width = 12, height = 6, file=paste0(rdaDir, "lac-reconstructed.pdf"))
par(oma=c(0,1.5,0,0))
layout(rbind(c(1:5), c(6:10), c(11,11,11,11,11)), heights = c(8,8,1))
for (ii in 1:10) {
  par(mar = c(4, 2.5, 1.75, 0.1))
  # plot(c(min(tvecsolve), max(tvecsolve)), c(min(c(recon[,i+1], xtrue[xtrue$time >=1,i+1], xsim.obs[,i+1])), max(c(recon[,i+1], xtrue[xtrue$time>=1,i+1], xsim.obs[,i+1]))), type = 'n', ylab='', xlab=paste('Component', i))
  # points(xsim[,1], xsim[,i+1], col='gray50')
  # lines(xtrue$time[xtrue$time >=1], xtrue[xtrue$time >=1,i+1], col="red")
  # lines(tvecsolve, recon[,i+1])
  # 
  # ourEst <- apply(oursPostExpX[,i,], 1, function(x) quantile(x, 0.5))
  # ourUB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.025)
  # ourLB <- apply(oursPostExpX[,i,], 1, quantile, probs = 0.975)
  times <- recon[,1]

  #plot(times, ourMed[,i], type="n", xlab="time", ylab=compnames[i]) #, ylim=c(ylim_lower[i], ylim_upper[i]))
  #if (ii ==7)
  #  plot(c(min(times), max(times)), c(min(c(ourLB[,ii], xtrue[,ii+1])), max(c(ourUB[,ii], xtrue[,ii+1]))), type = 'n', ylab='', xlab='', ylim=c(88,108))
  #else
  plot(c(min(times), max(times)), yaxisLims[,ii], type = 'n', ylab='', xlab='')
  
  mtext(compnames[ii], cex=1.25)
  polygon(c(times, rev(times)), c(ourUB[,ii], rev(ourLB[,ii])),
          col = "grey75", border = NA)
  lines(xtrue$time, xtrue[,ii+1], col="red", lwd=2)
  lines(times, ourMed[,ii], col="forestgreen", lwd=1.5)

  if (ii == 1)
    mtext("       Concentration (arb. unit)",side=2,line=0,outer=TRUE)
  if (ii == 8)
    title(xlab='Time (sec)', cex.lab = 1.5)
    
}
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)

legend("center",  c("truth", "median of all reconstructed trajectories", "95% interval from the 2.5 and 97.5 percentile of all reconstructed trajectories"),
       lty = c(1, 1, 0), lwd = c(2, 2, 0), bty = "n",
       col = c("red", "forestgreen", NA), fill = c(0, 0, "grey75"), text.width=c(0.05, 0.2, 0.25),
       border = c(0, 0, "grey75"), pch = c(NA, NA, 15), horiz = TRUE, cex = 1.25)


dev.off()


sink(paste0(rdaDir, "/result", expt, ".txt"))
   
library(xtable)

# theta posterior mean table
oursPostTheta <- sapply(oursPostTheta, identity, simplify = "array")

mean_est <- rbind(
  rowMeans(oursPostTheta[-1,1,])
)

sd_est <- rbind(
  apply(oursPostTheta[-1,1,], 1, sd)
)

rmse_est <- rbind(
  apply(oursPostTheta[-1,1,] - pram.true$theta[-1], 1, function(x) sqrt(mean(x^2)))
)


digit_precision = 4
printr <- function(x, precision=digit_precision) format(round(x, precision), nsmall=precision)
tablizeEstErr <- function(est, err){
  paste(format(round(est, digit_precision), nsmall=digit_precision), "\\pm", format(round(err, digit_precision), nsmall=digit_precision))
}

tab <- cbind(
  paste0("$k", 1:16, "$"), tablizeEstErr(mean_est[1,],sd_est[1,]), signif(rmse_est[1,],3)
)
tab <- data.frame(tab)
#rownames(tab) <- c("Method", "a", "b", "c", "d")
#colnames(tab) <- NULL
library(xtable)
print(xtable(tab), include.rownames=FALSE)

sink()
