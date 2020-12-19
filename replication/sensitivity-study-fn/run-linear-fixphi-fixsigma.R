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

outDir <- paste0("../results/linear-fill", filllevel, "-nobs", nobs_keep, "/")
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

temperature <- "heating"
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

plotPostSamplesFlexLinear(
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
sigma1 <- rbind(cbind(gpcov$Cinv, 0), 0)
sigma2 <- t(cbind(gpcov$mphi, -1)) %*% gpcov$Kinv %*% cbind(gpcov$mphi, -1)
dd <- diag(c(is.finite(xsim[,2]) * config$noise^2, 0))
ztilde <- c(xsim[,2], 0)
ztilde[is.na(ztilde)] <- 0

posterior_mean <- solve(sigma1 + sigma2 + dd, dd%*%ztilde)
posterior_variance <- solve(sigma1 + sigma2 + dd)

x0theta_id <- c(1, length(posterior_mean))
x0theta_mean <- posterior_mean[x0theta_id]
x0theta_variance <- posterior_variance[x0theta_id, x0theta_id]

library(mvtnorm)
plot(gpode$xsampled[id_plot,1,], gpode$theta[id_plot], xlab="x0", ylab="theta")
points(x0theta_mean[1], x0theta_mean[2], pch=20, cex=2, col=2)

plot.ellipse <- function (A, mu, r, n.points = 1000) {
  theta <- seq(0, 2 * pi, length = n.points)
  v <- rbind(r * cos(theta), r * sin(theta))
  ## transform for points on ellipse
  z <- backsolve(chol(A), v) + mu
  ## plot points
  lines(t(z), type = "l", col=2, lwd=2)
}

plot.ellipse(solve(x0theta_variance), x0theta_mean, -2*log(0.05))
plot.ellipse(solve(x0theta_variance), x0theta_mean, -2*log(0.32))
legend("topright", c("95% contour", "68% contour"), lty=1, col=2, lwd=2)
