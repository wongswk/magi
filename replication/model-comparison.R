# data A ----
load("../results/Michaelis-Menten-Va/summary-Michaelis-Menten-Va-fill0.5-noise0.04-phi132.1-datava-csv-time0to40-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_a <- new.env())
model_a <- as.list(model_a)
load("../results/Michaelis-Menten-Vb4p/summary-Michaelis-Menten-Vb4p-fill0.5-noise0.02-phi201.7-datava-csv-time0to40-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_b <- new.env())
model_b <- as.list(model_b)

noise_data <- 0.02

oos_rmse <- list()
oos_ppp <- list()
for(it in 1:100){
  xsim.obs <- read.csv(paste0("../results/Michaelis-Menten-Va/va_xsim_obs_seed",it,".csv"), row.names = 1)
  time_oos <- tail(xsim.obs$time, 3)
  idx_oos <- time_oos * 2 + 1
  xsim_oos <- tail(xsim.obs, 3)[,-1]
  xinfer_oos_a <- model_a$ours[[it]]$xpostmean[idx_oos, c(2,3)]
  err_oos_a <- xinfer_oos_a - xsim_oos
  
  xinfer_oos_b <- model_b$ours[[it]]$xpostmean[idx_oos, c(2,4)]
  err_oos_b <- xinfer_oos_b - xsim_oos
  
  err_mat <- rbind(sqrt(colMeans(err_oos_a^2)), sqrt(colMeans(err_oos_b^2)))
  rownames(err_mat) <- c("A", "B")
  oos_rmse[[it]] <- err_mat
  
  # Posterior Predictive P-value
  oos_ppp[[it]] <- sapply(c("model_a", "model_b"), function(model_name){
    oos_samples <- get(model_name)$outStorage[[it]]$oos_samples
    oos_samples <- oos_samples + rnorm(length(oos_samples), mean = 0, sd = noise_data)
    ppp <- rowMeans(apply(oos_samples, 1, function(x) x > xsim_oos))
    ppp <- pmin(ppp, 1-ppp)
    ppp <- matrix(ppp*2, ncol=2)
    ppp
  }, simplify = "array")
}

oos_rmse <- sapply(oos_rmse, identity, simplify = "array")
oos_ppp <- sapply(oos_ppp, identity, simplify = "array")

c1 <- rgb(0,0,1,1/4)
c2 <- rgb(1,0,0,1/4)
pdf("../results/histogram data A OOS RMSE PPP model A model B.pdf")
for (component in c("S", "P")){
  hgA <- hist(oos_rmse["A",component,], plot=FALSE)
  hgB <- hist(oos_rmse["B",component,], plot=FALSE)

  plot(hgA, col = c1, xlim = range(oos_rmse[,component,]), ylim = c(0,30), 
       main=paste0("OOS RMSE, data A, component ", component))
  plot(hgB, add = TRUE, col = c2)
  legend("topright", c("A", "B"), col=c(c1, c2), lty=1, lwd=20)
  hist(oos_rmse["A",component,] - oos_rmse["B",component,], main=paste0("OOS RMSE A - B, data B, component ", component))
  abline(v=0, col=2, lwd=3)
  
  tab <- rbind(summary(oos_rmse["A",component,]),
               summary(oos_rmse["B",component,]),
               summary(oos_rmse["A",component,] - oos_rmse["B",component,]))
  rownames(tab) <- c(paste0("data A - OOS RMSE of model A on component ", component),
                     paste0("data A - OOS RMSE of model B on component ", component),
                     paste0("data A - OOS RMSE difference of model A - B on component ", component))
  print(tab)
}

layout(matrix(1:6, ncol=2))
dimnames(oos_ppp)[[2]] <- c("S", "P")
for (component in c("S", "P")){
  for (it in 1:3){
    hgA <- hist(oos_ppp[it,component,"model_a",], plot=FALSE)
    hgB <- hist(oos_ppp[it,component,"model_b",], plot=FALSE)
    
    plot(hgA, col = c1, xlim = range(oos_ppp[it,component,,]), ylim = c(0,20), 
         main=paste0("OOS PPP, data A, component ", component, ", time ", tail(xsim.obs$time, 3)[it]))
    plot(hgB, add = TRUE, col = c2)
  }
}

dev.off()


# data B ----
load("../results/Michaelis-Menten-Va/summary-Michaelis-Menten-Va-fill0.5-noise0.04-phi132.1-datavb-csv-time0to40-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_a <- new.env())
model_a <- as.list(model_a)
load("../results/Michaelis-Menten-Vb4p/summary-Michaelis-Menten-Vb4p-fill0.5-noise0.02-phi201.7-datavb-csv-time0to40-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda", model_b <- new.env())
model_b <- as.list(model_b)
noise_data <- 0.01

oos_rmse <- list()
is_rmse <- list()
oos_ppp <- list()
xdesolveTRUE <- read.csv("../results/Michaelis-Menten-Vb4p.csv", row.names = 1)
obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
for(it in 1:100){
  xsim.obs <- read.csv(paste0("../results/Michaelis-Menten-Vb4p/vb_xsim_obs_seed",it,".csv"), row.names = 1)
  xsim.obs <- xsim.obs[obs_keep,]
  xsim.obs <- rbind(c(0, 1, 0), xsim.obs)
  time_oos <- tail(xsim.obs$time, 3)
  idx_oos <- time_oos * 2 + 1
  
  idx_is <- sapply(xsim.obs[1:(nrow(xsim.obs)-3),"time"], function(x) which(abs(x-xdesolveTRUE[,1]) < 1e-6))
  xdesolveTRUE.obs <- xdesolveTRUE[idx_is, c("S", "P")]
  idx_is <- sapply(xsim.obs[1:(nrow(xsim.obs)-3),"time"], function(x) which(abs(x-seq(0,70,0.5)) < 1e-6))
  xpostmean_a.obs <- model_a$ours[[it]]$xpostmean[idx_is, c(2,3)]
  xpostmean_b.obs <- model_b$ours[[it]]$xpostmean[idx_is, c(2,4)]
  is_err_mat <- rbind(
    sqrt(colMeans((xdesolveTRUE.obs - xpostmean_a.obs)^2)),
    sqrt(colMeans((xdesolveTRUE.obs - xpostmean_b.obs)^2))
  )
  
  rownames(is_err_mat) <- c("A", "B")
  is_rmse[[it]] <- is_err_mat
  
  xsim_oos <- tail(xsim.obs, 3)[,-1]
  xinfer_oos_a <- model_a$ours[[it]]$xpostmean[idx_oos, c(2,3)]
  err_oos_a <- xinfer_oos_a - xsim_oos
  
  xinfer_oos_b <- model_b$ours[[it]]$xpostmean[idx_oos, c(2,4)]
  err_oos_b <- xinfer_oos_b - xsim_oos
  
  err_mat <- rbind(sqrt(colMeans(err_oos_a^2)), sqrt(colMeans(err_oos_b^2)))
  rownames(err_mat) <- c("A", "B")
  oos_rmse[[it]] <- err_mat
  
  # Posterior Predictive P-value
  oos_ppp[[it]] <- sapply(c("model_a", "model_b"), function(model_name){
    oos_samples <- get(model_name)$outStorage[[it]]$oos_samples
    oos_samples <- oos_samples + rnorm(length(oos_samples), mean = 0, sd = noise_data)
    ppp <- rowMeans(apply(oos_samples, 1, function(x) x > xsim_oos))
    ppp <- pmin(ppp, 1-ppp)
    ppp <- matrix(ppp*2, ncol=2)
    ppp
  }, simplify = "array")
}

oos_rmse <- sapply(oos_rmse, identity, simplify = "array")
is_rmse <- sapply(is_rmse, identity, simplify = "array")
oos_ppp <- sapply(oos_ppp, identity, simplify = "array")

c1 <- rgb(0,0,1,1/4)
c2 <- rgb(1,0,0,1/4)
pdf("../results/histogram data B OOS RMSE PPP model A model B.pdf")
for (component in c("S", "P")){
  hgA <- hist(oos_rmse["A",component,], plot=FALSE)
  hgB <- hist(oos_rmse["B",component,], plot=FALSE)
  
  plot(hgA, col = c1, xlim = range(oos_rmse[,component,]), ylim = c(0,30), 
       main=paste0("OOS RMSE, data B, component ", component))
  plot(hgB, add = TRUE, col = c2)
  legend("topright", c("A", "B"), col=c(c1, c2), lty=1, lwd=20)
  hist(oos_rmse["A",component,] - oos_rmse["B",component,], main=paste0("OOS RMSE A - B, data B, component ", component))
  abline(v=0, col=2, lwd=3)
  
  tab <- rbind(summary(oos_rmse["A",component,]),
        summary(oos_rmse["B",component,]),
        summary(oos_rmse["A",component,] - oos_rmse["B",component,]))
  rownames(tab) <- c(paste0("data B - OOS RMSE of model A on component ", component),
                     paste0("data B - OOS RMSE of model B on component ", component),
                     paste0("data B - OOS RMSE difference of model A - B on component ", component))
  print(tab)
}

layout(matrix(1:6, ncol=2))
dimnames(oos_ppp)[[2]] <- c("S", "P")
for (component in c("S", "P")){
  for (it in 1:3){
    hgA <- hist(oos_ppp[it,component,"model_a",], plot=FALSE)
    hgB <- hist(oos_ppp[it,component,"model_b",], plot=FALSE)
    
    plot(hgA, col = c1, xlim = range(oos_ppp[it,component,,]), ylim = c(0,20), 
         main=paste0("OOS PPP, data A, component ", component, ", time ", tail(xsim.obs$time, 3)[it]))
    plot(hgB, add = TRUE, col = c2)
  }
}
avg_ppp <- apply(oos_ppp, 3:4, mean)

hgA <- hist(avg_ppp["model_a",], plot=FALSE)
hgB <- hist(avg_ppp["model_b",], plot=FALSE)
layout(1)

plot(hgA, col = c1, xlim = range(oos_ppp[it,component,,]), ylim = c(0,80), 
     main=paste0("OOS PPP, data A, component ", component, ", time ", tail(xsim.obs$time, 3)[it]))
plot(hgB, add = TRUE, col = c2)

hist(avg_ppp["model_a",] - avg_ppp["model_b",])
abline(v=0, col=2, lwd=3)

dev.off()

apply(is_rmse, 1:2, mean)


