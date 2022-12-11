library(magi)

outDir <- "../results/Michaelis-Menten-Inhibitor/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0(outDir, "hydrolysis.csv"))
matplot(realdata$t, realdata[,-1], type="b")

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = nrow(realdata),
    noise = c(NA, 0.02, 0.02, NA, NA),
    kernel = "generalMatern",
    seed = 12345,
    bandsize = 20,
    hmcSteps = 100,
    n.iter = 20001,
    stepSizeFactor = 0.01,
    linfillspace = 0.5, 
    t.end = 70,
    modelName = "Michaelis-Menten-Inhibitor"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(0.9, 0.75, 2.54, 10, 0.5),
  x0 = c(0.1, 1, 0, 0.05, 0),
  phi = cbind(c(0.1, 70), c(1, 30), c(1, 30), c(1, 30), c(1, 30)),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenInhibitorODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(3,4)], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = round(realdata$t / config$linfillspace) * config$linfillspace)
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))
xtest <- xsim

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
xsim.obs <- rbind(c(0, pram.true$x0), xsim.obs)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

## Linearly interpolate using fixed interval widths
fillC <- seq(0, config$t.end, by = config$linfillspace)
xsim <- data.frame(time = fillC)
xsim <- cbind(xsim, matrix(NaN, nrow = length(fillC), ncol = ncol(xsim.obs)-1 ))
for (i in 1:length(fillC)) {
  loc <- match( fillC[i], xsim.obs[, "time"])
  if (!is.na(loc))
    xsim[i,2:ncol(xsim)] = xsim.obs[loc,2:ncol(xsim)];
}

# cpp inference ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenInhibitorODE,
  fOdeDx=magi:::MichaelisMentenInhibitorDx,
  fOdeDtheta=magi:::MichaelisMentenInhibitorDtheta,
  thetaLowerBound=c(0,-100,0,0,-100),
  thetaUpperBound=c(Inf,Inf,Inf,Inf,Inf),  
  name="Michaelis-Menten-Inhibitor"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
xsim[1,5] <- NaN
for(j in c(1,2,3,5)){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0
  }
}


OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control = 
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, useFixedSigma=TRUE, verbose=TRUE))

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

gpode$lglik <- gpode$lp
pram.true$sigma <- sigma_fixed

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[2], ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
tail(gpode$theta)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)




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


stepSizeFactor <- rep(0.01, (nrow(xsim))*(length(pram.true$x0)-2) + length(dynamicalModelListReduced$thetaLowerBound) + length(pram.true$x0) - 2)
for(j in 1:3){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0
  }
}

reduced_result <- magi::MagiSolver(xsim[,c(2,3,4)], dynamicalModelListReduced, xsim$time, control = 
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  burninRatio = 0.5, phi = pram.true$phi[,1:3], sigma=sigma_fixed[1:3], useFixedSigma=TRUE, verbose=TRUE))

gpode <- reduced_result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, dynamicalModelListReduced$fOde(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))


