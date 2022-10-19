library(magi)

outDir <- "../results/Michaelis-Menten-Log1x/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0("../results/Michaelis-Menten/hydrolysis.csv"))
matplot(realdata$t, realdata[,-1], type="b")

# set up configuration if not already exist ------------------------------------
# if(!exists("config")){
config <- list(
  nobs = nrow(realdata),
  noise = c(0.01, 0.01, 0.01, 0.01),  # noise = c(0.01, 0.01, 0.01, 0.01), for fully observed case
  kernel = "generalMatern",
  seed = 123,
  bandsize = 100,
  hmcSteps = 100,
  n.iter = 20001,
  stepSizeFactor = 0.01,
  linfillspace = 0.5, 
  t.end = 70,
  modelName = "Michaelis-Menten-Log(1+x)"
)
# }


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(0.686, 0.01, 2.80),
  x0 = log(1+c(0.1, 1, 0, 0)),
  phi = cbind(c(1, 30), c(1, 60), c(1, 30), c(1, 10)),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenlog1xModelODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(3,5)], type="l", lty=1)

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
  fOde=magi:::MichaelisMentenlog1xModelODE,
  fOdeDx=magi:::MichaelisMentenlog1xModelDx,
  fOdeDtheta=magi:::MichaelisMentenlog1xModelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="Michaelis-Menten-Log(1+x)"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-6

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

xInitExogenous <- data.matrix(xsim[,-1])
for (j in c(2,4)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}
xInitExogenous[, 1] <- 0.1
xInitExogenous[-1, 3] <- 0.05

xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
# xInitExogenous <- NULL



OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control = 
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = config$stepSizeFactor,
                                  xInit = xInitExogenous, burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, discardBurnin=TRUE, useFixedSigma=TRUE,
                                  skipMissingComponentOptimization=TRUE, useMean=FALSE, useBand=TRUE))

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
gpode$theta <- cbind(gpode$theta, (gpode$theta[,2]+gpode$theta[,3])/gpode$theta[,1])
pram.true$theta <- c(pram.true$theta, (pram.true$theta[2]+pram.true$theta[3])/pram.true$theta[1])

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-fill", sum(config$linfillspace),"-noise", sum(config$noise, na.rm = TRUE), "-phi", sum(pram.true$phi), ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
tail(gpode$theta)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)
