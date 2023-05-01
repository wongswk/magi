# Replication script for Figure 9 (Effect of hyperparameter values phi)

library(magi)

seeds <- 100

# Run 100 repetitions with different random seeds. Can be done in parallel.
for (i in 1:seeds) {
  system(paste0("Rscript HIV-random-phi.R ", i))
}

# Read summary data and produce Figure 9
HIVsummary <- list()
for (seed in 1:seeds) {
  load(paste0("random-phi-", seed, ".rda"))
  #print(summary(HIVresult, sigma = TRUE))
  #print(HIVresult$phi)
  
  HIVsummary[[seed]] <- summary(HIVresult, digits = 6)
}

# 95% CI limits from original analysis
ci_lower <- c(34.1, 0.0982, 0.494, 943, 2.84, 2.79, 2.95, 2.23)
ci_upper <- c(37.7, 0.1150, 0.504, 973, 2.92, 3.71, 3.96, 37.80)

pdf("HIV-sens.pdf", width = 12, height = 4)
par(mfrow=c(1,5))
HIVpar.names <- c("lambda", "rho", "delta", "N", "c")
for (i in 1:length(HIVpar.names)) {
  par(mar = c(1, 2, 1.5, 0.5))
  estCI <- sapply(HIVsummary, function(x) x[,i])
  plot(1:seeds, xlim = c(0, seeds+1), ylim = c(0.975*min(estCI[2, ]), 1.025*max(estCI[3, ])),
       xaxt = 'n', xlab = '', ylab = '', type = 'n')
  segments(1:seeds, y0 = estCI[2, ], y1 = estCI[3, ], col = 1:seeds)
  mtext(HIVpar.names[i])
  points(1:seeds, estCI[1, ], col = 1:seeds, cex = 1)
  
  abline(h = ci_lower[i], lwd = 2, lty = 2)
  abline(h = ci_upper[i], lwd = 2, lty = 2)
}
dev.off()
