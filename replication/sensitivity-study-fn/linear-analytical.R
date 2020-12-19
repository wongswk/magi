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

temperature <- "notempering"
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

# comparing with analytical result

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

plot.ellipse <- function (A, mu, r, n.points = 1000, add=TRUE) {
  theta <- seq(0, 2 * pi, length = n.points)
  v <- rbind(r * cos(theta), r * sin(theta))
  ## transform for points on ellipse
  z <- backsolve(chol(A), v) + mu
  ## plot points
  if(add){
    lines(t(z), type = "l", col=2, lwd=2)  
  }else{
    plot(t(z), type = "l", col=2, lwd=2, xlab="x0", ylab="theta")
  }
}

samples_iid <- rmvnorm(2001, x0theta_mean, x0theta_variance)

# pdf(paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "analytical-contour.pdf"), height = 8, width = 8)
plot(samples_iid, pch=2, col=1, xlab="x0", ylab="theta")
points(x0theta_mean[1], x0theta_mean[2], pch=20, cex=2, col=2)
plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.95,df=2)))
plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.68,df=2)))
legend("topright", c("95% contour", "68% contour"), lty=1, col=2, lwd=2)
mtext(paste0(capture.output(round(x0theta_mean, 3), round(x0theta_variance, 5)), collapse = "\n"))
# dev.off()

plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.95,df=2)), add=FALSE)
points(x0theta_mean[1], x0theta_mean[2], pch=20, cex=2, col=2)
plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.68,df=2)))
legend("topright", c("95% contour", "68% contour"), lty=1, col=2, lwd=2)
mtext(paste0(capture.output(round(x0theta_mean, 3), round(x0theta_variance, 5)), collapse = "\n"))


prior_sd_x0 <- 2
prior_sd_theta <- 0.2

get_posterior <- function(xsim, temperature="notempering"){
  if(temperature == "heating"){
    config$priorTemperature <- nrow(xsim) / sum(is.finite(xsim[,2]))
    config$priorTemperatureObs <- 1
  }else if(temperature == "cooling"){
    config$priorTemperature <- 1
    config$priorTemperatureObs <- sum(is.finite(xsim[,2])) / nrow(xsim)
  }else if(temperature == "notempering"){
    config$priorTemperature <- 1
    config$priorTemperatureObs <- 1
  }else{
    stop("temperature must be 'heating' or 'cooling' or 'notempering'")
  }
  
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
  
  plot(NA, xlim = pram.true$x0 + c(-prior_sd_x0, prior_sd_x0), ylim=pram.true$theta + c(-prior_sd_theta, prior_sd_theta), xlab="x0", ylab="theta")
  plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.95,df=2)))
  points(x0theta_mean[1], x0theta_mean[2], pch=20, cex=2, col=2)
  plot.ellipse(solve(x0theta_variance), x0theta_mean, sqrt(qchisq(0.68,df=2)))
  legend("topright", c("95% contour", "68% contour"), lty=1, col=2, lwd=2)
  mtext(paste0(capture.output(round(x0theta_mean, 3), round(x0theta_variance, 5)), collapse = "\n"))
  return(list(x0theta_mean=x0theta_mean, x0theta_variance=x0theta_variance))
}

outDir <- paste0("../results/linear-nobs", nobs_keep, "/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
for(temperature in c("notempering", "heating", "cooling")){
  pdf(paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1],  "-", temperature, "-analytical-convergence.pdf"), height = 8, width = 8)
  for(each_fill_level in 0:6){
    get_posterior(insertNaN(xsim.obs, each_fill_level), temperature)
    legend("topleft", paste0("ndis =", nrow(insertNaN(xsim.obs, each_fill_level))))
  }
  dev.off()
}
