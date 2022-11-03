library(magi)

# try default phi
# what if we observe only after the drop (after t = 2)
# 

outDir <- "../results/Michaelis-Menten/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0(outDir, "hydrolysis.csv"))

if(!exists("externalFlag")){
  externalFlag <- FALSE
}

# set up configuration if not already exist ------------------------------------
if(!externalFlag){
  config <- list(
    nobs = nrow(realdata),
    noise = c(NaN, 0.02, 0.02),  # noise = c(0.01, 0.01, 0.01, 0.01), for fully observed case
    kernel = "generalMatern",
    seed = 123,
    bandsize = 40,
    hmcSteps = 100,
    n.iter = 5001,
    linfillspace = c(0.1, 0.2, 0.5), 
    linfillcut = c(2, 5), 
    t.end = 70,
    t.start = 0,
    t.truncate = 70,
    useMean = TRUE,
    phi = cbind(c(0.1, 70), c(1, 30), c(1, 30)),
    modelName = "Michaelis-Menten-Reduced"
  )
}

if(is.null(config$skip_visualization)){
  matplot(realdata$t, realdata[,-1], type="b")
}

# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list(
  theta=c(0.9, 0.75, 2.54),
  # x0 = c(0.08277011, 0.7500533, 0.2327168),
  x0 = c(0.1, 1, 0),
  phi = config$phi,
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenReducedODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
if(is.null(config$skip_visualization)){
  matplot(xtrue[, "time"], xtrue[, c(3,4)], type="l", lty=1)  
}


xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

if(length(config$linfillcut) == 0){
  xsim <- data.frame(time = round(realdata$t / config$linfillspace) * config$linfillspace)
}else{
  fill_seg <- c()
  startpoint = 0
  for(i in 1:length(config$linfillcut)){
    cutpoint <- config$linfillcut[i]
    fill_seg <- c(fill_seg, seq(startpoint, cutpoint, by = config$linfillspace[i]))
    startpoint <- cutpoint
  }
  fill_seg <- c(fill_seg, seq(startpoint, config$t.end, by = config$linfillspace[i+1]))
  xsim <- data.frame(time = fill_seg[sapply(realdata$t, function(t_each) which.min(abs(fill_seg - t_each)))])
}

xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))
xtest <- xsim

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
xsim.obs <- rbind(c(0, pram.true$x0), xsim.obs)

if(is.null(config$skip_visualization)){
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)
}

## Linearly interpolate using fixed interval widths
if(length(config$linfillcut) == 0){
  fillC <- seq(0, config$t.end, by = config$linfillspace)
}else{
  fillC <- fill_seg
}

xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}

xsim.obs <- xsim.obs[xsim.obs$time <= config$t.truncate, ]
xsim <- xsim[xsim$time <= config$t.truncate, ]
xsim.obs <- xsim.obs[xsim.obs$time >= config$t.start, ]
xsim <- xsim[xsim$time >= config$t.start, ]


# cpp inference ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenReducedODE,
  fOdeDx=magi:::MichaelisMentenReducedDx,
  fOdeDtheta=magi:::MichaelisMentenReducedDtheta,
  thetaLowerBound=c(0,-100,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="Michaelis-Menten-Reduced"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

# xInitExogenous <- data.matrix(xsim[,-1])
# for (j in c(2,4)){
#   xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
#   idx <- which(is.na(xInitExogenous[, j]))
#   xInitExogenous[idx, j] <- xInitExogenous[idx[1] - 1, j]
# }
# xInitExogenous[-1, 1] <- 0.1
# xInitExogenous[-1, 3] <- 0.05

xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
# xInitExogenous <- NULL

stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
# if(config$t.start == 0){
  for(j in 1:3){
    for(incre in 1:1){
      stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0  
    }
  }
# }
OursStartTime <- proc.time()[3]


result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control =
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  xInit = xInitExogenous, burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, discardBurnin=TRUE, useFixedSigma=TRUE,
                                  skipMissingComponentOptimization=TRUE, useMean=config$useMean, useBand=FALSE, priorTemperature=NULL))

OursTimeUsed <- proc.time()[3] - OursStartTime


gpode <- result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, dynamicalModelList$fOde(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = dynamicalModelList$fOde(pram.true$theta, data.matrix(xtrue[,-1]), xtrue$time)

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}
config$phi <- NULL

gpode$lglik <- gpode$lp
pram.true$sigma <- sigma_fixed
gpode$theta <- cbind(gpode$theta, (gpode$theta[,2]+gpode$theta[,3])/gpode$theta[,1])
pram.true$theta <- c(pram.true$theta, (pram.true$theta[2]+pram.true$theta[3])/pram.true$theta[1])

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-fill", sum(config$linfillspace),"-noise", sum(config$noise, na.rm = TRUE), "-phi", sum(pram.true$phi),"-useMean", config$useMean,"-time", config$t.start,"to", config$t.truncate, ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
tail(gpode$theta)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)

matplot(apply(gpode$xsampled[,,], 2:3, median), type="l")
