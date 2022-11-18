# data B ----
# run with less IS and more OOS
for(noise_scalar in c(0.02, 0.01, 0.005, 0.002, 0.001)){
  for(hold_out_size in seq(2, 14, 2)){
    
    obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
    obs_time = c(0.5, 1, 2.5, 3.5, 4.5, 5.5, 7, 8.5, 9.5, 11, 12, 13.5, 15, 
                 16, 18, 20, 21.5, 24, 27, 29.5, 32.5, 35.5, 39.5, 45, 55, 69)
    
    obs_time = obs_time[obs_keep]
    t.truncate = obs_time[length(obs_time) - hold_out_size]
    
    load(paste0("../results/Michaelis-Menten-Va/summary-Michaelis-Menten-Va-fill0.5-noise", 2 * noise_scalar,
                "-phi132.1-datavb-csv-time0to",t.truncate,"-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda"), model_a <- new.env())
    model_a <- as.list(model_a)
    load(paste0("../results/Michaelis-Menten-Vb4p/summary-Michaelis-Menten-Vb4p-fill0.5-noise", 2 * noise_scalar,
                "-phi201.7-datavb-csv-time0to",t.truncate,"-obs_keep3;5;7;9;10;..-linfillcut-time_changepoint0factor1.rda"), model_b <- new.env())
    model_b <- as.list(model_b)

    
    oos_rmse <- list()
    is_rmse <- list()
    oos_ppp <- list()
    xdesolveTRUE <- read.csv("../results/Michaelis-Menten-Vb4p.csv", row.names = 1)
    
    n_seed = min(length(model_a$ours), length(model_b$ours))
    for(it in 1:n_seed){
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
      
      # Posterior Predictive P-value
      oos_ppp[[it]] <- sapply(c("model_a", "model_b"), function(model_name){
        oos_samples <- get(model_name)$outStorage[[it]]$oos_samples
        oos_samples <- oos_samples + rnorm(length(oos_samples), mean = 0, sd = noise_scalar)
        ppp <- rowMeans(apply(oos_samples, 1, function(x) x > xsim_oos))
        avg_component <- apply(oos_samples, c(1,3), mean)
        ppp <- matrix(ppp, ncol=2)
        ppp <- rbind(ppp, rowMeans(apply(avg_component, 1, function(x) x > colMeans(xsim_oos))))
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
    pdf(paste0("../results/histogram data B OOS RMSE PPP model A model B noise ", noise_scalar, " hold_out ", hold_out_size, ".pdf"))
    for (component in c("S", "P")){
      hgA <- hist(oos_rmse["A",component,], plot=FALSE)
      hgB <- hist(oos_rmse["B",component,], plot=FALSE)
      
      plot(hgA, col = c1, xlim = range(oos_rmse[,component,]), ylim = c(0,30), 
           main=paste0("OOS RMSE, data B, component ", component))
      plot(hgB, add = TRUE, col = c2)
      legend("topright", c("A", "B"), col=c(c1, c2), lty=1, lwd=20)
      hist(oos_rmse["A",component,] - oos_rmse["B",component,], main=paste0("OOS RMSE A - B, data B, component ", component))
      abline(v=0, col=2, lwd=3)
      mtext(paste0("error rate = ", round(mean(oos_rmse["A",component,] - oos_rmse["B",component,] < 0)*100, 2), "%"))
      
      # tab <- rbind(summary(oos_rmse["A",component,]),
      #              summary(oos_rmse["B",component,]),
      #              summary(oos_rmse["A",component,] - oos_rmse["B",component,]))
      # rownames(tab) <- c(paste0("data B - OOS RMSE of model A on component ", component),
      #                    paste0("data B - OOS RMSE of model B on component ", component),
      #                    paste0("data B - OOS RMSE difference of model A - B on component ", component))
      # print(tab)
    }
    
    hgA <- hist(colSums(oos_rmse["A",,]), plot=FALSE)
    hgB <- hist(colSums(oos_rmse["B",,]), plot=FALSE)
    plot(hgA, col = c1, xlim = range(c(colSums(oos_rmse["B",,]), colSums(oos_rmse["A",,]))), ylim = c(0,max(c(hgA$counts, hgB$counts))), 
         main=paste0("OOS RMSE, data B, component both"))
    plot(hgB, add = TRUE, col = c2)
    legend("topright", c("A", "B"), col=c(c1, c2), lty=1, lwd=20)
    hist( colSums(oos_rmse["A",,] - oos_rmse["B",,]), main=paste0("OOS RMSE A - B, data B, component both"))
    abline(v=0, col=2, lwd=3)
    mtext(paste0("error rate = ", round(mean(colSums(oos_rmse["A",,] - oos_rmse["B",,]) < 0)*100, 2), "%"))
    
    layout(matrix(1:2, ncol=2))
    dimnames(oos_ppp)[[2]] <- c("S", "P")
    for (component in c("S", "P")){
      for (it in 1:(hold_out_size+1)){
        hgA <- hist(oos_ppp[it,component,"model_a",], plot=FALSE)
        hgB <- hist(oos_ppp[it,component,"model_b",], plot=FALSE)
        
        plot(hgA, col = c1, xlim = range(oos_ppp[it,component,,]), ylim = c(0,max(c(hgA$counts, hgB$counts))), 
             main=paste0("OOS PPP, data A, component ", component, ", time ", tail(xsim.obs$time, hold_out_size)[it]))
        plot(hgB, add = TRUE, col = c2)
        legend("topright", paste0(c("A", "B"), "; ", round(c(mean(oos_ppp[it,component,"model_a",] < 0.05), mean(oos_ppp[it,component,"model_b",] < 0.05))*100, 2), "% < 0.05"),
               col=c(c1, c2), lty=1, lwd=20)
      }
    }
    
    # noise level 0.01
    # OOS size 
    
    avg_ppp <- apply(oos_ppp[-(hold_out_size+1),,,], 3:4, mean)
    
    hgA <- hist(avg_ppp["model_a",], plot=FALSE)
    hgB <- hist(avg_ppp["model_b",], plot=FALSE)
    layout(1)
    
    plot(hgA, col = c1, xlim = c(0, 1), ylim = c(0,max(c(hgA$counts, hgB$counts))), 
         main=paste0("average OOS PPP data B"))
    plot(hgB, add = TRUE, col = c2)
    
    legend("topright", paste0(c("A", "B"), "; ", round(c(mean(avg_ppp["model_a",] < 0.05), mean(avg_ppp["model_b",] < 0.05))*100, 2), "% < 0.05"),
           col=c(c1, c2), lty=1, lwd=20)
    
    print("is_rmse")
    is_rmse_avg <- apply(is_rmse, 1:2, mean)
    colnames(is_rmse_avg) <- paste0("IS-RMSE-", colnames(is_rmse_avg))
    
    infoTab <- as.data.frame(is_rmse_avg)
    infoPerRow <- 6
    npanel <- ceiling(ncol(infoTab)/infoPerRow)
    tbls <- lapply(1:npanel, function(i) {
      gridExtra::tableGrob(infoTab[, ((i - 1) * infoPerRow + 
                                        1):min(ncol(infoTab), i * infoPerRow)])
    })
    do.call(gridExtra::grid.arrange, c(tbls, nrow = length(tbls)))

    dev.off()
  }
}
