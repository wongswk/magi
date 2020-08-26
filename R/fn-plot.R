library(gpds)

load("../results-from-backup/compare-withRamsay.rda")

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
# sink(paste0(outDirWenk, "/result.txt"))
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
layout(t(1:2))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, ours[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB[,i], rev(ourLB[,i])),
          col = "grey80", border = NA)
  lines(times, ours[[1]]$xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourMed[,i], lwd=3)
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
  
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
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
  
}
dev.off()


desolveWenk <- sapply(Wenk, function(x) x$xdesolvePM[,-1], simplify = "array")
WenkLB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.025))
WenkMed <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.5))
WenkUB <- apply(desolveWenk, c(1,2), function(x) quantile(x, 0.975))

times <- Wenk[[1]]$xdesolveTRUE[,1]

pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotWenk.pdf"))
layout(t(1:2))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, Wenk[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(WenkUB[,i], rev(WenkLB[,i])),
          col = "grey80", border = NA)
  lines(times, Wenk[[1]]$xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, WenkMed[,i], lwd=3)
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
  
}

dev.off()

pdf(width = 20, height = 1.2, file=paste0(outDirWenk, "legendGrey.pdf"))
par(mar=rep(0,4))
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth", "median of all reconstructed trajectories", "95% interval from the 2.5 and 97.5 percentile of all reconstructed trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "black", NA), fill=c(0, 0, "grey80"), pch=c(NA, NA, 15), x.intersp=c(2.5,2.5,0),
       border=c(0, 0, "grey80"), bty = "n", cex=1.8)
dev.off()


desolveDondel <- sapply(Dondel, function(x) x$xdesolvePM[,-1], simplify = "array")
DondelLB <- apply(desolveDondel, c(1,2), function(x) quantile(x, 0.025))
DondelMed <- apply(desolveDondel, c(1,2), function(x) quantile(x, 0.5))
DondelUB <- apply(desolveDondel, c(1,2), function(x) quantile(x, 0.975))
times <- Dondel[[1]]$xdesolveTRUE[,1]
pdf(width = 20, height = 5, file=paste0(outDirWenk, "plotDondel.pdf"))
layout(t(1:2))
for (i in 1:(ncol(xsim)-1)) {
  plot(times, Dondel[[1]]$xdesolveTRUE[,1+i], type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(DondelUB[,i], rev(DondelLB[,i])),
          col = "grey80", border = NA)
  lines(times, Dondel[[1]]$xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, DondelMed[,i], lwd=3)
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
  
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
layout(t(1:2))
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
dev.off()

pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxOursHoriz.pdf"))

layout(t(1:2))
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
dev.off()



xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxWenk.pdf"))
layout(t(1:2))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(WenkPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(WenkPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(WenkPostX[,i,], 1, quantile, probs = 0.975)
  times <- xsim.obs$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(xdesolveTRUE[,1], xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourEst, col="forestgreen", lwd=3)
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
}
dev.off()
# sink()

xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = xsim$time, func = odemodel$modelODE, parms = pram.true$theta)
pdf(width = 20, height = 5, file=paste0(outDirWenk, "posteriorxDondel.pdf"))
layout(t(1:2))
for (i in 1:(ncol(xsim)-1)) {
  ourEst <- apply(DondelPostX[,i,], 1, quantile, probs = 0.5)
  ourUB <- apply(DondelPostX[,i,], 1, quantile, probs = 0.025)
  ourLB <- apply(DondelPostX[,i,], 1, quantile, probs = 0.975)
  times <- xsim.obs$time
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  lines(xdesolveTRUE[,1], xdesolveTRUE[,1+i], col="red", lwd=4)
  lines(times, ourEst, col="forestgreen", lwd=3)
  if(i==1){
    mtext("V", cex=2)
  }else{
    mtext("R", cex=2)
  }
}
dev.off()

pdf(width = 20, height = 1.2, file=paste0(outDirWenk, "legendBlue.pdf"))
par(mar=rep(0,4))
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
legend("center", c("truth", "median of all inferred trajectories", "95% interval from the 2.5 and 97.5 percentile of all inferred trajectories"), lty=c(1,1,0), lwd=c(4,3,0),
       col = c("red", "forestgreen", NA), fill=c(0, 0, "skyblue"), pch=c(NA, NA, 15), x.intersp=c(2.5,2.5,0),
       border=c(0, 0, "skyblue"), bty = "n", cex=1.8)
dev.off()


pdf(width = 20, height = 0.5, file=paste0(outDirWenk, "header.pdf"))
par(mar=rep(0,4))
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
text(0.5, 0.5, "Ours", cex=2, font=2)
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
text(0.5, 0.5, "Wenk", cex=2, font=2)
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
text(0.5, 0.5, "Dondelinger", cex=2, font=2)
plot(c(0,1), c(0,1) ,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)
text(0.5, 0.5, "Ramsay", cex=2, font=2)
dev.off()
