rdaDir <- "../results/Michaelis-Menten-Va/"

for(phi2 in c(70, 120)){
  for(scenario in 0:3){
    
    phi = cbind(c(0.1, phi2), c(1, 30), c(1, 30))
    
    linfillspace = c(0.5)
    linfillcut = NULL
    phi_change_time = 0
    time_acce_factor = 1
    noise = c(NA, 0.02, 0.02)
    obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
    
    if(scenario == 0){
      obs_source = "va-csv"
      t.truncate = 70
    }else if (scenario == 1){
      obs_source = "vb-csv"
      t.truncate = 70
    }else if(scenario == 2){
      obs_source = "va-csv"
      t.truncate = 40
    }else if (scenario == 3){
      obs_source = "vb-csv"
      t.truncate = 40
    }
    
    config <- list(
      nobs = 26,
      noise = noise,
      kernel = "generalMatern",
      seed = 1,
      bandsize = 40,
      hmcSteps = 100,
      n.iter = 8001,
      linfillspace = linfillspace, 
      linfillcut = linfillcut,
      t.end = 70,
      t.start = 0,
      obs_start_time = 0,
      phi_change_time = phi_change_time,
      time_acce_factor = time_acce_factor,
      t.truncate = t.truncate,
      obs_keep = obs_keep,
      useMean = TRUE,
      phi = phi,
      skip_visualization = TRUE,
      obs_source = obs_source,
      modelName = "Michaelis-Menten-Va"
    )
    
    if(!is.null(config$linfillcut)){
      config$linfillcut <- paste(round(config$linfillcut, 2), collapse = ";")
      config$linfillspace <- paste(round(config$linfillspace, 2), collapse = ";")
    }else{
      config$linfillcut <- NULL
    }
    
    config$obs_keep <- paste(c(config$obs_keep[1:5], ".."), collapse = ";")
    {
    rm(list=setdiff(ls(), c("rdaDir", "subdirs", "config", "scenario", "phi2")))
    
    filename = paste0(rdaDir, "summary-", config$modelName,"-fill", config$linfillspace,"-noise", 
                      sum(config$noise, na.rm = TRUE), "-phi", sum(config$phi),
                      "-data", config$obs_source,
                      "-time", config$t.start,"to", config$t.truncate,
                      "-obs_keep", config$obs_keep, "-linfillcut", config$linfillcut,
                      "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor)
    single_rda_filename <- paste0(rdaDir, config$modelName,"-1-fill", config$linfillspace,"-noise",
                                  sum(config$noise, na.rm = TRUE), "-phi", sum(config$phi),
                                  "-data", config$obs_source,
                                  "-time", config$t.start,"to", config$t.truncate,
                                  "-obs_keep", config$obs_keep,
                                  "-linfillcut", config$linfillcut,
                                  "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
                                  ".rda")
    
    load(paste0(filename,".rda"))
    load(single_rda_filename)}
    if(config$obs_source == "va-csv"){
      xdesolveTRUE <- read.csv("../results/Michaelis-Menten-Va.csv", row.names = 1)
    }else{
      xdesolveTRUE <- read.csv("../results/Michaelis-Menten-Vb4p.csv", row.names = 1)
      xdesolveTRUE <- xdesolveTRUE[,c("time", "E", "S", "P")]
    }
    
    # load(paste0(rdaDir, rda_files[1]), envir = .GlobalEnv)
    rdaDir <- rdaDirSummary
    
    library(magi)
    library(xtable)
    
    # Average the posterior mean RMSEs for the different seeds
    starttime <- xsim$time[which(is.finite(xsim[,2]))[1]]
    xsim.obs <- xsim.obs[xsim.obs$time >= starttime,]
    rowId <- sapply(xsim.obs$time, function(x) which(abs(x-times) < 1e-6))
    
    for (i in 1:length(ours)) {
      rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xdesolveTRUE[,1]) < 1e-6))
      xdesolveTRUE.obs <- xdesolveTRUE[rowId,-1]
      rowId <- sapply(xsim.obs$time, function(x) which(abs(x-ours[[i]]$xdesolvePM[,1]) < 1e-6))
      xdesolvePM.obs <- ours[[i]]$xdesolvePM[rowId,-1]
      rowId <- sapply(xsim.obs$time, function(x) which(abs(x-xsim$time) < 1e-6))
      xpostmean.obs <- ours[[i]]$xpostmean[rowId,]
      
      ours[[i]]$rmseOdeExpPM <- sqrt(apply(((xdesolvePM.obs) - (xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
      ours[[i]]$rmseInferredTrajectory <- sqrt(apply(((xpostmean.obs) - (xdesolveTRUE.obs))^2, 2, mean, na.rm=TRUE))   # compared to true traj
    }
    rmse_reconstructed_orig <- round(apply(sapply(ours, function(x) x$rmseOdeExpPM), 1, mean), digits=4)
    print(rmse_reconstructed_orig)
    rmse_inferred_orig <- round(apply(sapply(ours, function(x) x$rmseInferredTrajectory), 1, mean), digits=4)
    print(rmse_inferred_orig)
    
    
    oursPostX <- sapply(oursPostX, identity, simplify = "array")
    
    ylim_lower <- rep(0, 10)
    ylim_upper <- apply(xdesolveTRUE, 2, max)[-1]
    
    pdf(width = 20, height = 5, file=paste0(filename, ".pdf"))
    # layout(rbind(c(1,2,3,4), c(5,5,5,5)), heights = c(5,1))
    
    matplot(xtrue[, "time"], (xtrue[, -1]), type="l", lty=1, col=c(1:6), xlab="time", ylab=NA)
    matplot(xsim.obs$time, (xsim.obs[,-1]), type="p", col=c(1:6), pch=20, add = TRUE)
    mtext('sample observations', cex=1.5)
    title("Michaelis Menten")
    # legend("topleft", c("true mRNA lacI", "true mRNA tetR", "true mRNA cI", "true protein lacI", "true protein tetR", "true protein cI",
    #                     "observed mRNA lacI", "observed mRNA tetR", "observed mRNA cI"), 
    #        lty=c(1,1,1,1,1,1,NA,NA,NA), pch=c(NA,NA,NA,NA,NA,NA,20,20,20), col=c(1:6, 1:3), cex=1.5)
    
    phiVisualization <- phiExogenous <- pram.true$phi
    compnames <- c("E", "S", "ES", "P")
    
    # smooth visualization with illustration
    
    id <- seq(1, nrow(xdesolveTRUE), by=1)
    xdesolveTRUE <- xdesolveTRUE[id,]
    
    xdesolvePM <- list()
    for (f in 1:length(rda_files)) {
      xdesolvePM[[f]] <- ours[[f]]$xdesolvePM
    }
    xdesolvePM_time <- (xdesolvePM[[1]][,1])
    xdesolvePM <- sapply(xdesolvePM, function(x) x[,-1], simplify = "array")
    
    id <- seq(1, nrow(xdesolvePM), by=2)
    ourXdesolveLB <- apply(xdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.025))
    ourXdesolveMed <- apply(xdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.5))
    ourXdesolveUB <- apply(xdesolvePM[id,,], c(1,2), function(x) quantile(x, 0.975))
    timesLUB <- xdesolvePM_time[id]
    
    for (i in 1:(ncol(xsim)-1)) {
      ourEst <- apply(oursPostX[,i,], 1, quantile, probs = 0.5)
      ourUB <- apply(oursPostX[,i,], 1, quantile, probs = 0.025)
      ourLB <- apply(oursPostX[,i,], 1, quantile, probs = 0.975)
      
      if(i != 1){
        ourEst <- (magi:::getMeanCurve(xsim$time, (ourEst), xdesolveTRUE[,1],
                                       t(phiVisualization[,i]), 0,
                                       kerneltype=config$kernel, deriv = FALSE))
        ourUB <- (magi:::getMeanCurve(xsim$time, (ourUB), xdesolveTRUE[,1],
                                      t(phiVisualization[,i]), 0,
                                      kerneltype=config$kernel, deriv = FALSE))
        ourLB <- (magi:::getMeanCurve(xsim$time, (ourLB), xdesolveTRUE[,1],
                                      t(phiVisualization[,i]), 0,
                                      kerneltype=config$kernel, deriv = FALSE))
      }else{
        ourEst <- approx(xsim$time, (ourEst), xdesolveTRUE[,1])$y
        ourUB <- approx(xsim$time, (ourUB), xdesolveTRUE[,1])$y
        ourLB <- approx(xsim$time, (ourLB), xdesolveTRUE[,1])$y
      }
      
      times <- xdesolveTRUE[,1]
      
      plot(times, ourEst, type="n", xlab="time", ylab=compnames[i], ylim=c(ylim_lower[i], ylim_upper[i]))
      
      polygon(c(times, rev(times)), c(ourUB, rev(ourLB)),
              col = "skyblue", border = "skyblue", lty = 1, density = 10, angle = -45)
      
      polygon(c(timesLUB, rev(timesLUB)), c(ourXdesolveUB[,i], rev(ourXdesolveLB[,i])),
              col = "grey80", border = "grey80", lty = 1, density = 10, angle = 45)
      
      lines(times, xdesolveTRUE[,1+i], col="red", lwd=4)
      lines(timesLUB, ourXdesolveMed[,i], lwd=3)
      lines(times, ourEst, col="forestgreen", lwd=3)
      mtext(compnames[i])
      # FIXME: impossible to have inferred trajectory and reconstructed trajectory differ at t=0
      # > mean(oursPostX[1,1,])
      # [1] 0.2286637
      # > mean(xdesolvePM[1,1,])
      # [1] 0.2286637
      
      # it is the start time bug
    }
    par(mar=rep(0,4))
    plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = TRUE)
    legend("center", c("truth", "median inferred trajectory", "median reconstructed trajectory", 
                       "CI on inferred trajectory", "CI on reconstructed trajectory"), lty=c(1,1,1,0,0), lwd=c(4,3,3,0,0),
           col = c("red", "forestgreen", "black", NA, NA), density=c(NA, NA, NA, 40, 40), fill=c(0, 0, 0, "skyblue", "grey80"), 
           border=c(0, 0, 0, "skyblue", "grey80"), angle=c(NA,NA,NA,-45,45), x.intersp=c(2.5,2.5,2.5,0, 0),  bty = "n", cex=1.8)
    dev.off()
    
    
    oursPostTheta <- sapply(oursPostTheta, identity, simplify = "array")
    
    printr <- function(x) format(round(x, 4), nsmall=4)
    tablizeEstErr <- function(est, err){
      paste(format(round(est, 4), nsmall=4), "\\pm", format(round(err, 4), nsmall=4))
    }
    
    mean_est <- rbind(
      rowMeans(oursPostTheta[,1,])
    )
    
    sd_est <- rbind(
      apply(oursPostTheta[,1,], 1, sd)
    )
    
    tab <- rbind(
      c("Ours", tablizeEstErr(mean_est[1,],sd_est[1,]))
    )
    tab <- data.frame(tab)
    colnames(tab) <- c("Method", c("k1", "K_{-1}", "k2", "km"))
    rownames(tab) <- NULL
    coverage <- rbind(
      printr(rowMeans((oursPostTheta[,2,] <= pram.true$theta) & (pram.true$theta <= oursPostTheta[,3,])))
    )
    coverage <- cbind(c("Ours"), coverage)
    colnames(coverage) <- c("Method", c("k1", "K_{-1}", "k2", "km"))
    
    tab <- cbind(c("truth", pram.true$theta), t(tab), t(coverage))
    
    theta_rmse <- sqrt(rowMeans((oursPostTheta[,1,] - pram.true$theta)^2))
    tab <- cbind(tab, c("theta rmse", printr(theta_rmse)))
    
    sink(paste0(filename, ".txt"))
    print("length(ours)=")
    print(length(ours))
    print("rmse_inferred_orig")
    print(rmse_inferred_orig)
    print("rmse_reconstructed_orig")
    print(rmse_reconstructed_orig)
    print(tab)
    # print(xtable(tab))
    sink()
  }
}
