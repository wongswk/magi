library(magi)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 15,
    noise = rep(0.001, 5), # 0.001 = low noise, 0.01 = high noise
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    bandsize = 40,
    hmcSteps = 100,
    n.iter = 20001,
    stepSizeFactor = 0.01,
    linfillspace = 0.5,  # discretization interval width (instead of fill level, since unevenly spaced observations)
    t.end = 100,
    modelName = "PTrans"
  )
}

config$ndis <- config$t.end / config$linfillspace + 1

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,100,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::ptransmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
#matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
#matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

## Linearly interpolate using fixed interval widths
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}

exoxInit <- sapply(2:ncol(xsim.obs), function(j)
  approx(xsim.obs[, "time"], xsim.obs[, j], xsim[, "time"])$y)


# cpp inference ----------------------------
ptransmodel <- list(
  name= config$modelName,
  fOde=magi:::ptransmodelODE,
  fOdeDx=magi:::ptransmodelDx,
  fOdeDtheta=magi:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(4,6)
)

outDir <- "../results/ptrans/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
config$priorTemperature <- config$ndis / config$nobs
config$priorTemperatureObs <- 1

OursStartTime <- proc.time()[3]

### Optimize phi first using equally spaced intervals of 1, i.e., 0,1...,100.
result <- magi::MagiSolver(exoxInit[xsim$time %in% 0:100,], ptransmodel, 0:100, control = list(
  nstepsHmc=config$hmcSteps, niterHmc=2, stepSizeFactor = config$stepSizeFactor, bandsize=config$bandsize, priorTemperature=config$priorTemperature))
sigmaUsed <- as.vector(result$sigma)

## NEW (Aug 11) ----- plug in sigma estimate and re-estimate phi
result <- magi::MagiSolver(exoxInit[xsim$time %in% 0:100,], ptransmodel, 0:100, control = list(
  nstepsHmc=config$hmcSteps, niterHmc=2, stepSizeFactor = config$stepSizeFactor, bandsize=config$bandsize, priorTemperature=config$priorTemperature, sigma=sigmaUsed))
phiUsed <- result$phi
sigmaUsed <- as.vector(result$sigma)

## final samples
result <- magi::MagiSolver(xsim[,-1], ptransmodel, xsim$time, control = list(
  nstepsHmc=config$hmcSteps, niterHmc=config$n.iter, stepSizeFactor = config$stepSizeFactor, bandsize=config$bandsize, priorTemperature=config$priorTemperature, sigma=sigmaUsed, phi=phiUsed, xInit=exoxInit))


OursTimeUsed <- proc.time()[3] - OursStartTime

gpode <- result
gpode$lglik <- gpode$lp
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, magi:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = magi:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, phiUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1],"-fill", config$linfillspace, ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-PTrans_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-PTrans_inferred_theta.csv"))

