
noise_scalar = 0.02
hold_out_size = 10
mtext_cex = 1.5
oos_bg_col = "lightpink"

obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
obs_time = c(0.5, 1, 2.5, 3.5, 4.5, 5.5, 7, 8.5, 9.5, 11, 12, 13.5, 15, 
             16, 18, 20, 21.5, 24, 27, 29.5, 32.5, 35.5, 39.5, 45, 55, 69)

obs_time = obs_time[obs_keep]
t.truncate = obs_time[length(obs_time) - hold_out_size]

load(paste0("../results/summary-Michaelis-Menten-Va-fill0.5-noise", 2 * noise_scalar,
            "-phi132.1-datavb-csv-time0to",t.truncate,"-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda"), model_a <- new.env())
model_a <- as.list(model_a)
load(paste0("../results/summary-Michaelis-Menten-Vb4p-fill0.5-noise", 2 * noise_scalar,
            "-phi201.7-datavb-csv-time0to",t.truncate,"-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda"), model_b <- new.env())
model_b <- as.list(model_b)


load("../results/Michaelis-Menten-Va/Michaelis-Menten-Va-1-fill0.5-noise0.04-phi132.1-datavb-csv-time0to20-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_1a <- new.env())
model_1a <- as.list(model_1a)
# FIXME model A and B phi should be consistent, either c(0.5, 30) or c(1, 30)
load("../results/Michaelis-Menten-Vb4p/Michaelis-Menten-Vb4p-1-fill0.5-noise0.04-phi201.7-datavb-csv-time0to20-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_1b <- new.env())
model_1b <- as.list(model_1b)


oos_rmse <- list()
is_rmse <- list()
oos_ppp <- list()
xdesolveTRUE <- read.csv("../results/Michaelis-Menten-Vb4p.csv", row.names = 1)

n_seed = min(length(model_a$ours), length(model_b$ours))
if(n_seed < 100){
  warning(paste("noise_scalar",noise_scalar,"hold_out_size",hold_out_size,"n_seed < 100"))
}

it = 1

xsim.obs <- read.csv(paste0("../results/Michaelis-Menten-Vb4p/noise",noise_scalar,"/vb_xsim_obs_seed",it,".csv"), row.names = 1)
xsim.obs <- xsim.obs[obs_keep,]
xsim.obs <- rbind(c(0, 1, 0), xsim.obs)
time_oos <- tail(xsim.obs$time, 3)

idx_all <- c(1, 6, 10, 15, 20, 23, 28, 31, 33, 37, 41, 44, 49, 55, 60, 66, 72, 80, 91, 111, 139)
idx_oos <- tail(idx_all, hold_out_size)  # only works for fill spacing = 0.5

idx_is <- sapply(xsim.obs[1:(nrow(xsim.obs)-hold_out_size),"time"], function(x) which(abs(x-xdesolveTRUE[,1]) < 1e-6))
xdesolveTRUE.obs <- xdesolveTRUE[idx_is, c("S", "P")]
idx_is <- sapply(xsim.obs[1:(nrow(xsim.obs)-hold_out_size),"time"], function(x) which(abs(x-seq(0,70,0.5)) < 1e-6))
xpostmean_a.obs <- model_a$ours[[it]]$xpostmean[idx_is, c(2,3)]
xpostmean_b.obs <- model_b$ours[[it]]$xpostmean[idx_is, c(2,4)]
is_err_mat <- rbind(
  sqrt(colMeans((xdesolveTRUE.obs - xpostmean_a.obs)^2)),
  sqrt(colMeans((xdesolveTRUE.obs - xpostmean_b.obs)^2))
)

rownames(is_err_mat) <- c("A", "B")
is_rmse[[it]] <- is_err_mat

xsim_oos <- tail(xsim.obs, hold_out_size)[,-1]
xinfer_oos_a <- model_a$ours[[it]]$xpostmean[idx_oos, c(2,3)]
err_oos_a <- xinfer_oos_a - xsim_oos

xinfer_oos_b <- model_b$ours[[it]]$xpostmean[idx_oos, c(2,4)]
err_oos_b <- xinfer_oos_b - xsim_oos

err_mat <- rbind(sqrt(colMeans(err_oos_a^2)), sqrt(colMeans(err_oos_b^2)))
rownames(err_mat) <- c("A", "B")
oos_rmse[[it]] <- err_mat


pdf(width = 15, height = 10, file="../results/MM-model-comparison.pdf")
# layout(cbind(c(7,1,1,6),c(2,2,4,4),c(3,3,5,5)))
layout(cbind(c(1,1,6,6),c(2,2,4,4),c(3,3,5,5)))
attach(model_1b)
matplot(xtrue[, "time"], (xtrue[, -1]), type="l", lty=1, col=0, xlab="time", ylab=NA)
matplot(xsim.obs$time[1:11], (xsim.obs[1:11,-1]), type="p", col=c(2,1), pch=19, add = TRUE)
matplot(xsim.obs$time[12:21], (xsim.obs[12:21,-1]), type="p", col=c(2,1), pch=5, add = TRUE)
mtext('observations', cex=2.5)

# legend("topright", c("observed [S] (in-sample)", "observed [P] (in-sample)",
#                      "observed [S] (hold-out)", "observed [P] (hold-out)"),
#        pch=c(19,19, 5, 5), col=c(2, 1, 2, 1), cex=1.5)

compnames = c("[E]", "[S]", "[ES]", "[P]")
ylim_lower <- rep(0, 10)
ylim_upper <- c(1,1,0.6,0.6)

for (i in c(4,2)) {
  phiVisualization <- pram.true$phi
  ourEst <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.5)
  ourUB <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.025)
  ourLB <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.975)
  
  
  ourEst <- (magi:::getMeanCurve(xsim$time, (ourEst), xdesolveTRUE[,1],
                                 t(phiVisualization[,i]), 0,
                                 kerneltype=config$kernel, deriv = FALSE))
  ourUB <- (magi:::getMeanCurve(xsim$time, (ourUB), xdesolveTRUE[,1],
                                t(phiVisualization[,i]), 0,
                                kerneltype=config$kernel, deriv = FALSE))
  ourLB <- (magi:::getMeanCurve(xsim$time, (ourLB), xdesolveTRUE[,1],
                                t(phiVisualization[,i]), 0,
                                kerneltype=config$kernel, deriv = FALSE))

  times <- xdesolveTRUE[,1]
  
  plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]), cex.lab=1.45)
  
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  
  idx_oos = (times > 20)
  polygon(c(times[idx_oos], rev(times[idx_oos])), 
          c(ourUB[idx_oos], rev(ourLB[idx_oos])),
          col = oos_bg_col, border = NA)
  
  lines(times, ourEst, col="forestgreen", lwd=3)
  
  if(i == 2){
    points(xsim.obs$time[1:11], (xsim.obs[1:11, 2]), type="p", col="red", pch=19)  
    points(xsim.obs$time[12:21], (xsim.obs[12:21, 2]), type="p", col="red", pch=5)
  }else{
    points(xsim.obs$time[1:11], (xsim.obs[1:11, 3]), type="p", col=1, pch=19)  
    points(xsim.obs$time[12:21], (xsim.obs[12:21, 3]), type="p", col=1, pch=5)
  }
  
  mtext(paste0("Component ",compnames[i], " inferred using model B"), cex=mtext_cex)
  if(i == 2){
    legend_loc = "topright"
  }else{
    legend_loc = "bottomright"
  }
  
  # FIXME: impossible to have inferred trajectory and reconstructed trajectory differ at t=0
  # > mean(oursPostX[1,1,])
  # [1] 0.2286637
  # > mean(xdesolvePM[1,1,])
  # [1] 0.2286637
  
  # it is the start time bug
}

detach(model_1b)
attach(model_1a)
for (i in c(3,2)) {
  phiVisualization <- pram.true$phi
  ourEst <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.5)
  ourUB <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.025)
  ourLB <- apply(gpode$xsampled[,,i], 2, quantile, probs = 0.975)
  
  
  ourEst <- (magi:::getMeanCurve(xsim$time, (ourEst), xdesolveTRUE[,1],
                                 t(phiVisualization[,i]), 0,
                                 kerneltype=config$kernel, deriv = FALSE))
  ourUB <- (magi:::getMeanCurve(xsim$time, (ourUB), xdesolveTRUE[,1],
                                t(phiVisualization[,i]), 0,
                                kerneltype=config$kernel, deriv = FALSE))
  ourLB <- (magi:::getMeanCurve(xsim$time, (ourLB), xdesolveTRUE[,1],
                                t(phiVisualization[,i]), 0,
                                kerneltype=config$kernel, deriv = FALSE))
  
  times <- xdesolveTRUE[,1]
  
  if(i == 2){
    legend_loc = "topright"
    comp = "[S]"
  }else{
    legend_loc = "bottomright"
    comp = "[P]"
  }
  
  plot(times, ourEst, type="n", xlab="time", ylab=comp, ylim=c(ylim_lower[i], ylim_upper[i]), cex.lab=1.45)
  
  polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)
  
  idx_oos = (times > 20)
  polygon(c(times[idx_oos], rev(times[idx_oos])), 
          c(ourUB[idx_oos], rev(ourLB[idx_oos])),
          col = oos_bg_col, border = NA)
  
  lines(times, ourEst, col="forestgreen", lwd=3)
  
  if(i == 2){
    points(xsim.obs$time[1:11], (xsim.obs[1:11, 2]), type="p", col="red", pch=19)  
    points(xsim.obs$time[12:21], (xsim.obs[12:21, 2]), type="p", col="red", pch=5)
  }else{
    points(xsim.obs$time[1:11], (xsim.obs[1:11, 3]), type="p", col=1, pch=19)  
    points(xsim.obs$time[12:21], (xsim.obs[12:21, 3]), type="p", col=1, pch=5)
  }
  
  mtext(paste0("Component ",comp, " inferred using model A"), cex=mtext_cex)

}

par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)

legend("center", c("observed [S] (in-sample)", "observed [P] (in-sample)",
                   "observed [S] (hold-out)", "observed [P] (hold-out)",
                   "posterior mean", "95% posterior interval for in-sample", "95% posterior interval for hold-out"), 
       lty=c(0,0,0,0,1,0,0), lwd=c(0,1,0,1,3,0,0),
       col = c(1,1,"red","red", "forestgreen", NA, NA), fill=c(0,0,0,0, 0,"skyblue",oos_bg_col),
       border=c(0,0,0,0, 0, "skyblue",oos_bg_col), pch=c(19,5,19,5, NA, 15, 15), cex=2.1)

dev.off()
