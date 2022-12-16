# Summarize the results
library(magi)


rdaDir <- "../results/repressilator-gene-regulation-log/51obs-fill1-noise0.3/"
pdf_files <- list.files(rdaDir)
rda_files <- pdf_files[grep("repressilator-gene-regulation-log.*\\.rda", pdf_files)]
pdf_files <- pdf_files[grep("repressilator-gene-regulation-log.*\\.pdf", pdf_files)]
  
  
config <- list()
config$modelName <- "repressilator-gene-regulation-log"
config$noise <- rep(0.3, 6)

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
  
  ourTime[jj] <- OursTimeUsed 
  
  xsampledexp <- exp(gpode$xsampled)
  oursPostExpX[[jj]] <- cbind(
    apply(xsampledexp, 2:3, mean),
    apply(xsampledexp, 2:3, function(x) quantile(x, 0.025)),
    apply(xsampledexp, 2:3, function(x) quantile(x, 0.975))
  )
  
  
}

print(paste0(rdaDir,"repressilator-log_summary.rda"))
save.image(file=paste0(rdaDir,"repressilator-log_summary.rda"))

xtrue <- xtrue[xtrue$time >= 1,]

jj <- 1 #sample(1:100, 1)
load(paste0(rdaDir, rda_files[jj]))

pdf(width = 12, height = 6, file=paste0(rdaDir, "rep-inferred-example.pdf"))
compnames <- c("m_lacI", "m_tetR", "m_cI", expression(paste("p_lacI (", bold("unobserved"), ")")), expression(paste("p_tetR (", bold("unobserved"), ")")), expression(paste("p_cI (", bold("unobserved"), ")")))
layout(rbind(c(1,2,3), c(4,5,6), c(7,7,7)), heights = c(8,8,1))
for (ii in 1:6) {
  
  par(mar = c(4, 4.5, 1.75, 0.1))
  ourEst <- oursPostExpX[[jj]][,ii]
  ourEst <- exp(magi:::getMeanCurve(xsim$time, log(ourEst), xtrue[,1],
                                    t(phiExogenous[,ii]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE))
  
  ourUB <- oursPostExpX[[jj]][,12+ii]
  ourUB <- exp(magi:::getMeanCurve(xsim$time, log(ourUB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))
  
  
  ourLB <- oursPostExpX[[jj]][,6+ii]
  ourLB <- exp(magi:::getMeanCurve(xsim$time, log(ourLB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))
  plot( c(min(xtrue$time),max(xtrue$time)), c(min(ourLB), min(max(ourUB),175)), type='n', xlab='', ylab='')
    
  polygon(c(xtrue[,1], rev(xtrue[,1])), c(ourUB, rev(ourLB)),
         col = "skyblue", border = NA)    
  # polygon(c(xsim[,1], rev(xsim[,1])), c(ourUB, rev(ourLB)),
  #        col = "skyblue", border = NA)    

  if (ii == 1)
    title(ylab='mRNA concentration', cex.lab = 1.5)

  if (ii == 4)
    title(ylab='Protein concentration', cex.lab = 1.5)

  if (ii == 5)
    title(xlab='Time (mRNA lifetimes)', cex.lab = 1.5)

  lines(xtrue[, "time"], exp(xtrue[,ii+1]),col='red', lwd=2)
  lines(xtrue[,1], ourEst, col='forestgreen', lwd=1.5)
  mtext(compnames[ii], cex=1.25)
  if (ii <= 3) points(xsim.obs$time[-1], exp(xsim.obs[-1,ii+1]), col='black', pch=16)
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


pdf(width = 12, height = 3.1, file=paste0(rdaDir, "rep-inferred-param-example.pdf"))
par.names <- c( expression(alpha[0]), expression(alpha), "n", expression(beta))
par(mfrow=c(1,4))
for (ii in 1:4) {
  if (ii == 1) par(oma=c(0,1.5,0,0))
  par(mar = c(2.5, 2.5, 2, 0.75))
  den <- density(gpode$theta[,ii])
  plot(den, main='', xlab = '', ylab = '', type='n')

  value1 <- quantile(gpode$theta[,ii], 0.025)
  value2 <- quantile(gpode$theta[,ii], 0.975)
  
  l <- min(which(den$x >= value1))
  h <- max(which(den$x < value2))
  
  polygon(c(den$x[c(l, l:h, h)]),
          c(0, den$y[l:h], 0),
          col = "grey75", border=NA)
  abline(v=pram.true$theta[ii], col='red', lwd =2)
      
  lines(den)
    
  
  if (ii == 1) mtext(text='Posterior density',side=2,line=0,outer=TRUE)
  mtext(par.names[ii], cex=1.25)
}
dev.off()



# Inferred trajectories
oursPostExpX <- simplify2array(oursPostExpX)
#xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = c(0,xsim$time), func = odemodel$modelODE, parms = pram.true$theta)
#xdesolveTRUE[,-1] <- exp(xdesolveTRUE[,-1])
pdf(width = 12, height = 6, file=paste0(rdaDir, "rep-posteriorxOurs-paper.pdf"))

layout(rbind(c(1,2,3), c(4,5,6), c(7,7,7)), heights = c(8,8,1))
for (ii in 1:6) {
  par(mar = c(4, 4.5, 1.75, 0.1))
  ourEst <- apply(oursPostExpX[,ii,], 1, function(x) quantile(x, 0.5))
  ourEst <- exp(magi:::getMeanCurve(xsim$time, log(ourEst), xtrue[,1],
                                    t(phiExogenous[,ii]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE))
    
  ourLB <- apply(oursPostExpX[,ii,], 1, quantile, probs = 0.025)
  ourLB <- exp(magi:::getMeanCurve(xsim$time, log(ourLB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))

  ourUB <- apply(oursPostExpX[,ii,], 1, quantile, probs = 0.975)
  ourUB <- exp(magi:::getMeanCurve(xsim$time, log(ourUB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))
  
  plot( c(min(xtrue$time),max(xtrue$time)), c(min(ourLB), min(max(ourUB),175)), type='n', xlab='', ylab='')  

  if (ii == 1)
    title(ylab='mRNA concentration', cex.lab = 1.5)
  
  if (ii == 4)
    title(ylab='Protein concentration', cex.lab = 1.5)
  
  if (ii == 5)
    title(xlab='Time (mRNA lifetimes)', cex.lab = 1.5)
  
  #plot(xtrue$time, ourEst, type="n", xlab="time", ylab=compnames[i]) #, ylim=c(ylim_lower[i], ylim_upper[i]))
  mtext(compnames[ii], cex=1.25)
  polygon(c(xtrue$time, rev(xtrue$time)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(xtrue$time, exp(xtrue[,1+ii]), col="red", lwd=2)
  lines(xtrue$time, ourEst, col="forestgreen", lwd=1.5)
}

par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
#legend("center", c("truth", "median of all inferred trajectories", "95% interval from the 2.5 and 97.5 percentile of all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
#       col = c("red", "forestgreen", NA), fill=c(0, 0,"skyblue"), text.width=c(0, 0.15, 0.01), bty = "n",
#       border=c(0, 0, "skyblue"), pch=c(NA, NA, 15), cex=1.2, horiz=TRUE)

legend("center",  c("truth", "median of all inferred trajectories", "95% interval from the 2.5 and 97.5 percentile of all inferred trajectories"),
       lty = c(1, 1, 0), lwd = c(2, 2, 0), bty = "n",
       col = c("red", "forestgreen", NA), fill = c(0, 0, "skyblue"), text.width=c(0.05, 0.2, 0.25),
       border = c(0, 0, "skyblue"), pch = c(NA, NA, 15), horiz = TRUE, cex = 1.25)
dev.off()

sink(paste0(rdaDir, "/result.txt"))
  
library(xtable)

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


digit_precision = 4
printr <- function(x, precision=digit_precision) format(round(x, precision), nsmall=precision)
tablizeEstErr <- function(est, err){
  paste(format(round(est, digit_precision), nsmall=digit_precision), "\\pm", format(round(err, digit_precision), nsmall=digit_precision))
}

tab <- cbind(
  1:4, tablizeEstErr(mean_est[1,],sd_est[1,]), rmse_est[1,]
)
tab <- data.frame(tab)
#rownames(tab) <- c("Method", "a", "b", "c", "d")
#colnames(tab) <- NULL
library(xtable)
print(xtable(tab), include.rownames=FALSE)

sink()
