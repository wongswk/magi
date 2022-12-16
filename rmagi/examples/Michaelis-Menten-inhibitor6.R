library(magi)

outDir <- "../results/Michaelis-Menten-Inhibitor/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

obs.times <- c(2.5, 4.5, 7, 9.5, 11, 13.5, 15, 16, 18, 20)
test.times <- c(21.5, 24, 27, 29.5, 32.5, 35.5, 39.5, 45, 55, 69)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = length(obs.times), #nrow(realdata),
    noise = c(NA, 0.02, 0.02, NA, NA, NA),
    kernel = "generalMatern",
    #seed = 12,
    seed = 669097609, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    #seed = 488500816,
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    linfillspace = 0.5, 
    t.end = 70,
    modelName = "MM-Inhibitor"
  )
}

# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(0.9, 0.75, 2.54, 1, 0.5),
  x0 = c(0.1, 1, 0, 0.2, 0, 0),
  phi = cbind(c(0.1, 70), c(1, 30), c(1, 30), c(1, 70), c(1, 70), c(1, 70)),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenInhibitor6ODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(3,4,5,6,7)], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = round(c(obs.times,test.times) / config$linfillspace) * config$linfillspace)
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

xtestDS <- xsim

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

# Divide into train/test
xtest <- xsim[xsim$time %in% test.times,]
xsim <- xsim[xsim$time %in% obs.times,]

xsim.obs <- rbind(c(0, pram.true$x0), xsim)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

## Linearly interpolate using fixed interval widths
xsim <- setDiscretization(rbind(xsim.obs, c(config$t.end, rep(NaN,ncol(xsim)-1))), by=config$linfillspace)
xsim[1,5] <- NaN  # do not observe initial I value

# cpp inference ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenInhibitor6ODE,
  fOdeDx=magi:::MichaelisMentenInhibitor6Dx,
  fOdeDtheta=magi:::MichaelisMentenInhibitor6Dtheta,
  thetaLowerBound=c(0,-100,0,0,-100),
  thetaUpperBound=c(Inf,Inf,Inf,Inf,Inf),  
  name="Michaelis-Menten-Inhibitor6"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtestDS[,-1]), pram.true$theta, xtestDS$time)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

# Fix initial values except for I
stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
for(j in c(1,2,3,5,6)){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0
  }
}


xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
xInitExogenous[,1] <- 0.1
xInitExogenous[,2] <- 1
xInitExogenous[,3] <- 0
xInitExogenous[,4] <- 0.1
xInitExogenous[,5] <- 0
xInitExogenous[,6] <- 0

OursStartTime <- proc.time()[3]

result2 <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control = 
                             list(xInit = xInitExogenous,bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  positiveSystem = TRUE, skipMissingComponentOptimization = TRUE, burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, useFixedSigma=TRUE, verbose=TRUE))

OursTimeUsed <- proc.time()[3] - OursStartTime

plot(result2)

gpode <- result2
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, dynamicalModelList$fOde(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))


xMean <- apply(gpode$xsampled, c(2, 3), mean)
fit_inhib <- xMean[gpode$tvec %in% xsim.obs[,"time"],2:3]
pred_inhib <- xMean[gpode$tvec %in% xtest[,"time"],2:3]
sse_train_inhib <- sum((fit_inhib - xsim.obs[,3:4])^2)
sse_test_inhib <-  sum((pred_inhib - xtest[,3:4])^2)

apply(gpode$xsampled[,1,], 2, median)
print(apply(gpode$theta, 2, median))

ourEst_inhib <- apply(gpode$xsampled, c(2, 3), mean)
ourLB_inhib <- apply(gpode$xsampled, c(2, 3), function(x) quantile(x, 0.025))
ourUB_inhib <- apply(gpode$xsampled, c(2, 3), function(x) quantile(x, 0.975))


#### compare to vanilla MM
# cpp inference ----------------------------
dynamicalModelListReduced <- list(
  fOde=magi:::MichaelisMentenReducedODE,
  fOdeDx=magi:::MichaelisMentenReducedDx,
  fOdeDtheta=magi:::MichaelisMentenReducedDtheta,
  thetaLowerBound=c(0,-100,0),
  thetaUpperBound=c(Inf,Inf,Inf),  
  name="Michaelis-Menten-Reduced"
)


stepSizeFactor <- rep(0.01, (nrow(xsim))*(length(pram.true$x0)-3) + length(dynamicalModelListReduced$thetaLowerBound) + length(pram.true$x0) - 3)
for(j in 1:3){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0
  }
}

reduced_result <- magi::MagiSolver(xsim[,c(2,3,4)], dynamicalModelListReduced, xsim$time, control = 
                             list(xInit = xInitExogenous[,1:3], bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  positiveSystem = TRUE, skipMissingComponentOptimization = TRUE, burninRatio = 0.5, phi = pram.true$phi[,1:3], sigma=sigma_fixed[1:3], useFixedSigma=TRUE, verbose=TRUE))

plot(reduced_result)

gpode <- reduced_result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, dynamicalModelListReduced$fOde(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))


xMean <- apply(gpode$xsampled, c(2, 3), mean)
fit_vanil <- xMean[gpode$tvec %in% xsim.obs[,"time"],2:3]
pred_vanil <- xMean[gpode$tvec %in% xtest[,"time"],2:3]
sse_train_vanil <- sum((fit_vanil - xsim.obs[,3:4])^2)
sse_test_vanil <-  sum((pred_vanil - xtest[,3:4])^2)

ourEst_vanil <- apply(gpode$xsampled, c(2, 3), mean)
ourLB_vanil <- apply(gpode$xsampled, c(2, 3), function(x) quantile(x, 0.025))
ourUB_vanil <- apply(gpode$xsampled, c(2, 3), function(x) quantile(x, 0.975))


pdf(paste0(outDir, config$modelName,"-",config$seed,"-nobs",config$nobs,"-noise", config$noise[2], "-prediction.pdf"), width=12, height=8)

# Visualization
compnames <- c("", "[S]", "[P]")

layout(cbind(c(1,1,6,6),c(2,2,4,4),c(3,3,5,5)))
par(mar = c(4, 4.5, 1.75, 0.1))

matplot(xtrue[, "time"], (xtrue[, -1]), type="n", lty=1, col=0, xlab='', ylab='mM')
title(xlab="Time (min)", line=2, cex.lab=1)
abline(v = max(obs.times), col="grey", lty=2, lwd=2)
matplot(xsim.obs$time, (xsim.obs[,3:4]), type="p", col=c(1,2), pch=19, add = TRUE)
matplot(xtest$time, xtest[,3:4], type="p", col=c(1,2), pch=5, add = TRUE)

mtext('observations', line = 0.3)

for (ii in 3:2) {
  
  par(mar = c(4, 4.5, 1.75, 0.1))
  ourEstp <- magi:::getMeanCurve(xsim$time, ourEst_inhib[,ii], xtrue[,1],
                                 t(pram.true$phi[,ii]), 0,
                                 kerneltype=config$kernel, deriv = FALSE)
  
  ourUBp <- magi:::getMeanCurve(xsim$time, ourUB_inhib[,ii], xtrue[,1],
                                t(pram.true$phi[,ii]), 0,
                                kerneltype=config$kernel, deriv = FALSE)
  
  
  ourLBp <- magi:::getMeanCurve(xsim$time, ourLB_inhib[,ii], xtrue[,1],
                                t(pram.true$phi[,ii]), 0,
                                kerneltype=config$kernel, deriv = FALSE)
  plot( c(min(xtrue$time),max(xtrue$time)), c(min(ourLBp), min(max(ourUBp),175)), type='n', xlab='', ylab='mM')
  title(xlab="Time (min)", line=2, cex.lab=1)
  abline(v = max(obs.times), col="grey", lty=2, lwd=2)
  
  polygon(c(xtrue[xtrue[,1] <= max(obs.times),1], rev(xtrue[xtrue[,1] <= max(obs.times),1])), c(ourUBp[xtrue[,1] <= max(obs.times)], rev(ourLBp[xtrue[,1] <= max(obs.times)])),
          col = "skyblue", border = NA)    
  polygon(c(xtrue[xtrue[,1] > max(obs.times),1], rev(xtrue[xtrue[,1] > max(obs.times),1])), c(ourUBp[xtrue[,1] > max(obs.times)], rev(ourLBp[xtrue[,1] > max(obs.times)])),
          col = "peachpuff", border = NA)      
  
  # lines(xtrue[, "time"], xtrue[,ii+1],col='red', lwd=2)
  lines(xtrue[,1], ourEstp, col='forestgreen', lwd=1.5)
  mtext(paste(compnames[ii], "inferred from inhibitor model"), line = 0.3)

  if(compnames[ii] == "[P]"){
    point_col = "red"
  }else{
    point_col = "black"
  }
  points(xsim$time, xsim[,ii+1], col=point_col, pch=16)
  points(xtest$time, xtest[,ii+1], col=point_col, pch=5)
}

for (ii in 3:2) {
  
  par(mar = c(4, 4.5, 1.75, 0.1))
  ourEstp <- magi:::getMeanCurve(xsim$time, ourEst_vanil[,ii], xtrue[,1],
                                 t(pram.true$phi[,ii]), 0,
                                 kerneltype=config$kernel, deriv = FALSE)
  
  ourUBp <- magi:::getMeanCurve(xsim$time, ourUB_vanil[,ii], xtrue[,1],
                                t(pram.true$phi[,ii]), 0,
                                kerneltype=config$kernel, deriv = FALSE)
  
  ourLBp <- magi:::getMeanCurve(xsim$time, ourLB_vanil[,ii], xtrue[,1],
                                t(pram.true$phi[,ii]), 0,
                                kerneltype=config$kernel, deriv = FALSE)
  plot( c(min(xtrue$time),max(xtrue$time)), c(min(ourLBp), min(max(ourUBp),175)), type='n', xlab='', ylab='mM')
  title(xlab="Time (min)", line=2, cex.lab=1)
  abline(v = max(obs.times), col="grey", lty=2, lwd=2)  
  
  polygon(c(xtrue[xtrue[,1] <= max(obs.times),1], rev(xtrue[xtrue[,1] <= max(obs.times),1])), c(ourUBp[xtrue[,1] <= max(obs.times)], rev(ourLBp[xtrue[,1] <= max(obs.times)])),
          col = "skyblue", border = NA)    
  polygon(c(xtrue[xtrue[,1] > max(obs.times),1], rev(xtrue[xtrue[,1] > max(obs.times),1])), c(ourUBp[xtrue[,1] > max(obs.times)], rev(ourLBp[xtrue[,1] > max(obs.times)])),
          col = "peachpuff", border = NA)    
  
    
  # lines(xtrue[, "time"], xtrue[,ii+1],col='red', lwd=2)
  lines(xtrue[,1], ourEstp, col='forestgreen', lwd=1.5)
  mtext(paste(compnames[ii], "inferred from M-M model"), line = 0.3)

  if(compnames[ii] == "[P]"){
    point_col = "red"
  }else{
    point_col = "black"
  }
  points(xsim$time, xsim[,ii+1], col=point_col, pch=16)
  points(xtest$time, xtest[,ii+1], col=point_col, pch=5)
  
}

# mtext(paste0("Inhibitor: SSE(train) = ", round(sse_train_inhib,3),
#              ", SSE(test) = ", round(sse_test_inhib,3),
#              "       M-M: SSE(train) = ", round(sse_train_vanil,3),
#              ", SSE(test) = ", round(sse_test_vanil,3))
#              ,side=1,line=1,outer=TRUE)

# par(mar = rep(0, 4))
# plot(1, type = 'n', xaxt = 'n', yaxt = 'n',
#      xlab = NA, ylab = NA, frame.plot = FALSE)
# 
# legend("top", c("truth", "inferred trajectory",
#                    "95% interval in training"),
#        lty = c(1, 1, 0), lwd = c(2, 2, 0), bty = "n",
#        col = c("red", "forestgreen", NA), fill = c(0, 0, "skyblue"),
#        border = c(0, 0, "skyblue"), pch = c(NA, NA, 15), horiz = TRUE, cex=1.1)
# legend("bottom", c("95% interval in prediction", "noisy observations for training", "noisy observations for prediction"),
#        lty = c(0, 0, 0), lwd = c(0, 1, 1), bty = "n",
#        col = c(NA, "black", "black"), fill = c("peachpuff", 0, 0),
#        border = c("peachpuff", 0, 0), pch = c(15, 16, 5), horiz = TRUE, cex=1.1)
par(mar=rep(0,4))
plot(1,type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, frame.plot = FALSE)

oos_bg_col = "peachpuff"

legend("center", c("observed noisy [S] for training", "observed noisy [S] for prediction",
                   "observed noisy [P] for training", "observed noisy [P] for prediction",
                   "inferred trajectory", "95% interval in training", "95% interval in prediction"), 
       lty=c(0,0,0,0,1,0,0), lwd=c(0,1,0,1,3,0,0),
       col = c(1,1,"red","red", "forestgreen", NA, NA), fill=c(0,0,0,0, 0,"skyblue",oos_bg_col),
       border=c(0,0,0,0, 0, "skyblue",oos_bg_col), pch=c(19,5,19,5, NA, 15, 15), cex=1.7)

dev.off()

cat(sse_train_inhib, sse_test_inhib, sse_train_vanil, sse_test_vanil, 
    file = paste0(outDir, config$modelName,"-",config$seed,"-nobs",config$nobs,"-noise", config$noise[2], "-sse.txt"))

sqrt(sse_train_inhib / nrow(xsim.obs)) # RMSE should be close to 0.02 if fit well
sqrt(sse_test_inhib / nrow(xtest)) # RMSE

sqrt(sse_train_vanil / nrow(xsim.obs)) # RMSE
sqrt(sse_test_vanil / nrow(xtest)) # RMSE
