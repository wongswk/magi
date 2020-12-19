library(magi)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0){
  filllevel <- as.numeric(args[1])
  seed <- scan("fn-seeds.txt")[as.numeric(args[2])]
  nobs_keep <- as.numeric(args[3])
}else{
  seed <- 117015794 # (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9
  filllevel <- 2
  nobs_keep <- 41
}

temperature <- "heating"
outDir <- paste0("../results/linear-fill", filllevel, "-nobs", nobs_keep, "-", temperature, "/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#' FIXME HMC in R has bug for random seed


# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = nobs_keep,
    noise = 1.5,
    kernel = "generalMatern",
    seed = seed,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = filllevel,
    t.end = 20,
    modelName = "linear",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = TRUE,
    max.epoch = 1
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperatureObs <- 1

if(config$loglikflag == "withmeanBand"){
  config$useMean = TRUE
  config$useBand = TRUE
}else if(config$loglikflag == "band"){
  config$useMean = FALSE
  config$useBand = TRUE
}else if(config$loglikflag == "withmean"){
  config$useMean = TRUE
  config$useBand = FALSE
}else if(config$loglikflag == "usual"){
  config$useMean = FALSE
  config$useBand = FALSE
}

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(1),
  x0 = c(-2),
  # phi=rbind(c(2.5, 0.75),
  #           c(1.5, 3.00)),
  phi=as.matrix(c(10, 3)),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=241)


fOde <- function(theta, x) {
  as.matrix(rep(theta, length(x)))
}

fOdeDx <- function (theta, x) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  resultDx
}

fOdeDtheta  <- function (theta, x) {
  resultDtheta <- array(1, c(nrow(x), length(theta), ncol(x)))
  resultDtheta
}
  
linearmodel <- list(
  fOde=fOde,
  fOdeDx=fOdeDx,
  fOdeDtheta=fOdeDtheta,
  thetaLowerBound=c(-Inf),
  thetaUpperBound=c(Inf),
  name="linear"
)

modelODE <- function(t, state, parameters) {
  list(as.vector(fOde(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,20,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- insertNaN(xsim.obs,config$filllevel)

# tempering must be disabled when comparing with analytical results
if(temperature == "heating"){
  config$priorTemperature <- config$ndis / config$nobs  
  config$priorTemperatureObs <- 1
}else if(temperature == "cooling"){
  config$priorTemperature <- 1
  config$priorTemperatureObs <- config$nobs / config$ndis
}else if(temperature == "notempering"){
  config$priorTemperature <- 1
  config$priorTemperatureObs <- 1
}else{
  stop("temperature must be 'heating' or 'cooling' or 'notempering'")
}

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

OursStartTime <- proc.time()[3]

samplesCpp <- magi:::solveMagiRcpp(
  yFull = data.matrix(xsim[,-1]),
  odeModel = linearmodel,
  tvecFull = xsim$time,
  sigmaExogenous = config$noise,
  phiExogenous = pram.true$phi,
  xInitExogenous = xInitExogenous,
  thetaInitExogenous = matrix(nrow=0,ncol=0),
  muExogenous = matrix(nrow=0,ncol=0),
  dotmuExogenous = matrix(nrow=0,ncol=0),
  priorTemperatureLevel = config$priorTemperature,
  priorTemperatureDeriv = config$priorTemperature,
  priorTemperatureObs = config$priorTemperatureObs,
  kernel = config$kernel,
  nstepsHmc = config$hmcSteps,
  burninRatioHmc = config$burninRatio,
  niterHmc = config$n.iter,
  stepSizeFactorHmc = config$stepSizeFactor,
  nEpoch = config$max.epoch,
  bandSize = config$bandsize,
  useFrequencyBasedPrior = config$useFrequencyBasedPrior,
  useBand = config$useBand,
  useMean = config$useMean,
  useScalerSigma = config$useScalerSigma,
  useFixedSigma = config$useFixedSigma,
  verbose = TRUE)

OursTimeUsed <- proc.time()[3] - OursStartTime

phiUsed <- samplesCpp$phi
samplesCpp <- samplesCpp$llikxthetasigmaSamples

samplesCpp <- samplesCpp[,,1]

out <- samplesCpp[-1,1,drop=FALSE]
xCpp <- matrix(out[1:length(data.matrix(xsim[,-1])), 1], ncol=ncol(xsim)-1)
thetaCpp <- out[(length(xCpp)+1):(length(xCpp) + length(linearmodel$thetaLowerBound)), 1]
sigmaCpp <- tail(out[, 1], ncol(xsim)-1)

matplot(xsim$time, xCpp, type="l", add=TRUE)

llikId <- 1
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(linearmodel$thetaLowerBound))
sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim)-1)


burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta=t(samplesCpp[thetaId, -(1:burnin), drop=FALSE]),
              xsampled=array(t(samplesCpp[xId, -(1:burnin)]),
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=samplesCpp[llikId,-(1:burnin)],
              sigma = t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, linearmodel$fOde(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = linearmodel$fOde(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(phiUsed[,j], 2), collapse = "; ")
}

plotPostSamplesFlex1d <- function (filename, xtrue, dotxtrue, xsim, gpode, param, config, odemodel = NULL) {
  npostplot <- config$npostplot
  if (is.null(npostplot)) {
    npostplot <- 20
  }
  xpostmean <- apply(gpode$xsampled, 2:3, mean)
  if (!is.null(odemodel)) {
    times <- sort(unique(round(c(odemodel$times, xsim$time, 
                                 xtrue[, "time"]), 7)))
    xdesolveTRUE <- deSolve::ode(y = param$x0, times = times, 
                                 func = odemodel$modelODE, parms = param$theta)
    mapId <- which.max(gpode$lglik)
    ttheta <- gpode$theta[mapId, ]
    tx0 <- gpode$xsampled[mapId, 1, ]
    xdesolveMAP <- deSolve::ode(y = tx0, times = times, 
                                func = odemodel$modelODE, parms = ttheta)
    ttheta <- colMeans(gpode$theta)
    tx0 <- mean(gpode$xsampled[, 1, ])
    xdesolvePM <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, 
                               parms = ttheta)
    rowId <- sapply(xsim$time, function(x) which(abs(x - 
                                                       times) < 1e-06))
    xdesolveTRUE.obs <- xdesolveTRUE[rowId, -1]
    xdesolveMAP.obs <- xdesolveMAP[rowId, -1]
    xdesolvePM.obs <- xdesolvePM[rowId, -1]
    xdesolveSamples <- parallel::mclapply(1:16, function(dummy) {
      mapId <- sample(1:length(gpode$lglik), 1)
      ttheta <- gpode$theta[mapId, ]
      tx0 <- gpode$xsampled[mapId, 1, ]
      deSolve::ode(y = tx0, times = odemodel$times, func = odemodel$modelODE, 
                   parms = ttheta)
    }, mc.cores = 8)
    rmseTrue <- sqrt(mean((xdesolveTRUE.obs - xsim[, -1])^2, na.rm = TRUE))
    rmseWholeGpode <- sqrt(mean((xpostmean - xsim[, -1])^2, na.rm = TRUE))
    rmseOdeMAP <- sqrt(mean((xdesolveMAP.obs - xsim[, -1])^2, na.rm = TRUE))
    rmseOdePM <- sqrt(mean((xdesolvePM.obs - xsim[, -1])^2, na.rm = TRUE))
    rmselist <- list(true = paste0(round(rmseTrue, 3), collapse = "; "), 
                     wholeGpode = paste0(round(rmseWholeGpode, 3), collapse = "; "), 
                     odeMAP = paste0(round(rmseOdeMAP, 3), collapse = "; "), 
                     odePM = paste0(round(rmseOdePM, 3), collapse = "; "))
    config$rmse <- rmselist
  }
  id.max <- which.max(gpode$lglik)
  config$noise <- paste(round(config$noise, 3), collapse = "; ")
  infoTab <- as.data.frame(config)
  pdf(filename, width = 8, height = 8)
  if (all(c("gridExtra", "gridBase") %in% rownames(installed.packages()))) {
    infoPerRow <- 6
    npanel <- ceiling(ncol(infoTab)/infoPerRow)
    tbls <- lapply(1:npanel, function(i) {
      gridExtra::tableGrob(infoTab[, ((i - 1) * infoPerRow + 
                                        1):min(ncol(infoTab), i * infoPerRow)])
    })
    do.call(gridExtra::grid.arrange, c(tbls, nrow = length(tbls)))
  }
  matplot(xtrue$time, xtrue[, -1], type = "l", lty = 1)
  matplot(xsim$time, xsim[, -1], type = "p", pch = 20, add = TRUE)
  id.plot <- seq(1, nrow(gpode$theta), length = npostplot)
  id.plot <- unique(as.integer(id.plot))
  id.plot <- unique(c(id.max, id.plot))
  for (j in 1:(ncol(xsim) - 1)) {
    matplot(xtrue$time, cbind(xtrue[, j + 1], dotxtrue[, 
                                                       j]), type = "l", lty = 1, col = c(2, 1), ylab = paste0("component-", 
                                                                                                              j), main = "full posterior")
    points(xsim$time, xsim[, j + 1], col = 2)
    matplot(xsim$time, t(gpode$xsampled[id.plot, , j]), 
            col = "skyblue", add = TRUE, type = "p", lty = 1, 
            pch = 20)
    matplot(xsim$time, t(gpode$xsampled[id.plot, , j]), 
            col = "skyblue", add = TRUE, type = "l", lty = 1)
    matplot(xsim$time, t(gpode$fode[id.plot, , j]), col = "grey", 
            add = TRUE, type = "p", lty = 1, pch = 20)
    matplot(xsim$time, t(gpode$fode[id.plot, , j]), col = "grey", 
            add = TRUE, type = "l", lty = 1)
    if (!is.null(odemodel) && !is.null(odemodel$curCov)) {
      lines(xsim$time, odemodel$curCov[[j]]$mu, col = "forestgreen", 
            lwd = 2)
      lines(xsim$time, odemodel$curCov[[j]]$dotmu, col = "darkgreen", 
            lwd = 2)
    }
    matplot(xtrue$time, cbind(xtrue[, j + 1], dotxtrue[, 
                                                       j]), type = "l", lty = 1, col = c(2, 1), ylab = paste0("component-", 
                                                                                                              j), main = "full posterior", add = TRUE)
  }
  if (!is.null(odemodel) && !is.null(odemodel$curCov)) {
    layout(1:2)
    for (j in 1:(ncol(xsim) - 1)) {
      matplot(xtrue$time, cbind(xtrue[, j + 1] - approx(xsim$time, 
                                                        odemodel$curCov[[j]]$mu, xtrue$time)$y, dotxtrue[, 
                                                                                                         j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, 
                                                                                                                     xtrue$time)$y), type = "l", lty = 1, col = c(2, 
                                                                                                                                                                  1), ylab = paste0("component-", j), main = "full posterior")
      points(xsim$time, xsim[, j + 1] - odemodel$curCov[[j]]$mu, 
             col = 2)
      matplot(xsim$time, t(gpode$xsampled[id.plot, , j]) - 
                odemodel$curCov[[j]]$mu, col = "skyblue", add = TRUE, 
              type = "p", lty = 1, pch = 20)
      matplot(xsim$time, t(gpode$xsampled[id.plot, , j]) - 
                odemodel$curCov[[j]]$mu, col = "skyblue", add = TRUE, 
              type = "l", lty = 1)
      matplot(xtrue$time, cbind(xtrue[, j + 1] - approx(xsim$time, 
                                                        odemodel$curCov[[j]]$mu, xtrue$time)$y, dotxtrue[, 
                                                                                                         j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, 
                                                                                                                     xtrue$time)$y), type = "l", lty = 1, col = c(2, 
                                                                                                                                                                  1), ylab = paste0("component-", j), main = "full posterior")
      matplot(xsim$time, t(gpode$fode[id.plot, , j]) - 
                odemodel$curCov[[j]]$dotmu, col = "grey", add = TRUE, 
              type = "p", lty = 1, pch = 20)
      matplot(xsim$time, t(gpode$fode[id.plot, , j]) - 
                odemodel$curCov[[j]]$dotmu, col = "grey", add = TRUE, 
              type = "l", lty = 1)
      matplot(xtrue$time, cbind(xtrue[, j + 1] - approx(xsim$time, 
                                                        odemodel$curCov[[j]]$mu, xtrue$time)$y, dotxtrue[, 
                                                                                                         j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, 
                                                                                                                     xtrue$time)$y), type = "l", lty = 1, col = c(2, 
                                                                                                                                                                  1), ylab = paste0("component-", j), main = "full posterior", 
              add = TRUE)
    }
  }
  layout(matrix(1:4, 2, byrow = TRUE))
  for (i in 1:length(param$theta)) {
    hist(gpode$theta[, i], main = letters[i])
    abline(v = param$theta[i], col = 2)
    abline(v = quantile(gpode$theta[, i], c(0.025, 0.975)), 
           col = 3)
  }
  if (length(gpode$sigma) < length(gpode$lglik)) {
    gpode$sigma <- matrix(gpode$sigma, nrow = 1)
  }
  for (sigmaIt in 1:ncol(gpode$sigma)) {
    hist(gpode$sigma[, sigmaIt], xlim = range(c(gpode$sigma[, 
                                                            sigmaIt], param$sigma[sigmaIt])), main = paste("sigma for component", 
                                                                                                           sigmaIt))
    abline(v = param$sigma, col = 2)
    abline(v = quantile(gpode$sigma[, sigmaIt], c(0.025, 
                                                  0.975)), col = 3)
  }
  layout(matrix(1:4, 2, byrow = TRUE))
  for (i in 1:length(param$theta)) {
    plot.ts(gpode$theta[, i], main = letters[i])
    abline(h = param$theta[i], col = 2)
    abline(h = quantile(gpode$theta[, i], c(0.025, 0.975)), 
           col = 3)
  }
  layout(1:2)
  plot(gpode$lglik, type = "l")
  hist(gpode$lglik)
  layout(1)
  xpostmedian <- apply(gpode$xsampled, 2:3, median)
  for (j in 1:(ncol(xsim) - 1)) {
    matplot(xtrue$time, cbind(xtrue[, j + 1], dotxtrue[, 
                                                       j]), type = "l", lty = 1, col = c(2, 1), ylab = paste0("component-", 
                                                                                                              j), main = "maximum a posterior")
    points(xsim$time, xsim[, j + 1], col = 2)
    points(xsim$time, gpode$xsampled[id.max, , j], col = "skyblue", 
           pch = 20)
    lines(xsim$time, gpode$xsampled[id.max, , j], col = "skyblue")
    points(xsim$time, gpode$fode[id.max, , j], col = "grey", 
           pch = 20)
    lines(xsim$time, gpode$fode[id.max, , j], col = "grey")
    lines(xsim$time, xpostmean[, j], type = "b", col = "green")
    lines(xsim$time, xpostmedian[, j], type = "b", col = "brown")
    legend("topleft", c("MAP", "Post-Mean", "Post-Median"), 
           col = c("skyblue", "green", "brown"), lty = 1)
  }
  if (!is.null(odemodel)) {
    layout(1)
    matplot(odemodel$xtrue[, "time"], odemodel$xtrue[, -1], 
            type = "l", lty = 1, add = FALSE, lwd = 3)
    for (xdesolveSampleEach in xdesolveSamples) {
      matplot(xdesolveSampleEach[, "time"], xdesolveSampleEach[, 
                                                               -1], type = "l", lty = 4, add = TRUE)
    }
    matplot(xdesolveMAP[, "time"], xdesolveMAP[, -1], type = "l", 
            lty = 3, add = TRUE, lwd = 2)
    matplot(xdesolvePM[, "time"], xdesolvePM[, -1], type = "l", 
            lty = 2, add = TRUE, lwd = 2)
    matplot(xsim$time, xsim[, -1], type = "p", col = 1:(ncol(xsim) - 
                                                          1), pch = 20, add = TRUE)
    legend("topleft", lty = c(1, 3, 2), legend = c("true", 
                                                   "MAP", "PosteriorMean"))
  }
  if (!is.null(odemodel) && !is.null(odemodel$curCov)) {
    curCov <- odemodel$curCov
    layout(1:(length(curCov) + 1))
    gpode$mxode <- sapply(1:length(curCov), function(j) curCov[[j]]$mphi %*% 
                            (t(gpode$xsampled[, , j]) - curCov[[j]]$mu) + curCov[[j]]$dotmu, 
                          simplify = "array")
    gpode$mxode <- aperm(gpode$mxode, c(2, 1, 3))
    gpode$odeErr <- gpode$fode - gpode$mxode
    postmeanOdeErr <- apply(gpode$odeErr, 2:3, mean)
    matplot(postmeanOdeErr, lty = 1, type = "l", main = "postmeanOdeErr")
    for (j in 1:length(curCov)) plot(xsim$time, cumsum(postmeanOdeErr[, 
                                                                      j]), lty = 1, type = "l", main = "cumsum postmeanOdeErr")
  }
  dev.off()
}


plotPostSamplesFlex1d(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_theta.csv"))


# comparing with analytical result
id_plot <- as.integer(seq(1, dim(gpode$xsampled)[1], length.out = 2001))
plot(gpode$theta[id_plot], gpode$xsampled[id_plot,1,], xlab = "theta", ylab="x0")

xtime <- xsim$time
gpcov <- calCov(pram.true$phi, 
                as.matrix(dist(xtime)),
                -sign(outer(xtime,xtime,'-')),
                kerneltype = "generalMatern",
                bandsize = config$bandsize)
sigma1 <- rbind(cbind(gpcov$Cinv / config$priorTemperature, 0), 0)
sigma2 <- t(cbind(gpcov$mphi, -1)) %*% (gpcov$Kinv / config$priorTemperature) %*% cbind(gpcov$mphi, -1)
dd <- diag(c(is.finite(xsim[,2]) * config$noise^(-2) / config$priorTemperatureObs, 0))

ztilde <- c(xsim[,2], 0)
ztilde[is.na(ztilde)] <- 0

posterior_mean <- solve(sigma1 + sigma2 + dd, dd%*%ztilde)
posterior_variance <- solve(sigma1 + sigma2 + dd)

x0theta_id <- c(1, length(posterior_mean))
x0theta_mean <- posterior_mean[x0theta_id]
x0theta_variance <- posterior_variance[x0theta_id, x0theta_id]

library(mvtnorm)

plot.ellipse <- function (A, mu, r, n.points = 1000) {
  theta <- seq(0, 2 * pi, length = n.points)
  v <- rbind(r * cos(theta), r * sin(theta))
  ## transform for points on ellipse
  z <- backsolve(chol(A), v) + mu
  ## plot points
  lines(t(z), type = "l", col=2, lwd=2)
}

samples_iid <- rmvnorm(length(id_plot), x0theta_mean, x0theta_variance)

pdf(paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "analytical-contour.pdf"), height = 8, width = 8)
plot(samples_iid, pch=2, col=2, xlab="x0", ylab="theta")
points(gpode$xsampled[id_plot,1,], gpode$theta[id_plot])
points(x0theta_mean[1], x0theta_mean[2], pch=20, cex=2, col=2)
plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.95,df=2)))
plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.68,df=2)))
legend("topright", c("95% contour", "68% contour"), lty=1, col=2, lwd=2)
mtext(paste0(capture.output(round(x0theta_mean, 3), round(x0theta_variance, 5)), collapse = "\n"))
dev.off()

save.image(paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))
