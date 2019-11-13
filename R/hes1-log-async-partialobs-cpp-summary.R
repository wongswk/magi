# Summarize the results
library(gpds)

pdf_files <- list.files("~/Workspace/DynamicSys/results/cpp/good/")
pdf_files <- pdf_files[grep("\\.pdf", pdf_files)]
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

save.image(file=paste0(rdaDir,"hes1log_summary.rda"))

# Average the posterior mean RMSEs for the different seeds
rmse.table <- round(apply(sapply(ours, function(x) x$rmseOdePM), 1, mean), digits=4)
print(rmse.table)

logscale <- FALSE
# Make the figures showing Ours using ODE solver results
# use the same axis limits for both methods for easier visual comparison
# on log scale
ylim_lower <- c(-1, -1, -2)
ylim_upper <- c(3, 2, 5)

# on original scale
if(!logscale){
  ylim_lower <- c(0, 0, 0)
  ylim_upper <- c(10, 4, 20)
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

pdf(width = 20, height = 5, file=paste0(rdaDir, "/plotOurs.pdf"))
par(mfrow=c(1, ncol(xsim)-1))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
  lines(times, ourMed[,i])
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
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, xdesolveTRUE[,1+i], col="red", lwd=1)
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

