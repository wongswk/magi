library(magi)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 33,
    noise = c(0.15,0.15,0.15),
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 2e4,
    stepSizeFactor = 0.01,
    filllevel = 0,
    modelName = "Hes1-log"
  )
}

# initialize global parameters, true x, simulated x ----------------------------
outDir <- "../results/hes1log/"
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = log(c(1.438575, 2.037488, 17.90385)),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::hes1logmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- xtrue

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}
xsim$X3 <- NaN
xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
xsim.obs$X1[seq(2,nrow(xsim.obs),by=2)] <- NaN
xsim.obs$X2[seq(1,nrow(xsim.obs),by=2)] <- NaN

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- setDiscretization(xsim.obs,config$filllevel)

# cpp inference ----------------------------
hes1logmodel <- list(
  name="Hes1-log",
  fOde=magi:::hes1logmodelODE,
  fOdeDx=magi:::hes1logmodelDx,
  fOdeDtheta=magi:::hes1logmodelDtheta,
  thetaLowerBound=rep(0,7),
  thetaUpperBound=rep(Inf,7)
)

# heating up with (3, 3, 1) temperature ------------------------------------
#' there are 99 sampled X's (3 x 33) and only 33 actual observations.  
#' Using the logic in our other two models, the default temperatures for Hes1 should be (3, 3, 1).

result <- magi::MagiSolver(xsim[,-1], hes1logmodel, xsim$time, control = list(nstepsHmc=config$hmcSteps, niterHmc=config$n.iter, stepSizeFactor = config$stepSizeFactor, sigma=pram.true$sigma, useFixedSigma=TRUE))
gpode <- result
gpode$lglik <- gpode$lp

gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, magi:::hes1logmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = magi:::hes1logmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)


for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}

save.image(paste0(outDir, config$modelName,"-",config$seed,"-heating.rda"))

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-heating.pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)


write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-hes1log_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-hes1log_inferred_theta.csv"))
