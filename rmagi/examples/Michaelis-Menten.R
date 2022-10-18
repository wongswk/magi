library(magi)

outDir <- "../results/Michaelis-Menten/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0(outDir, "hydrolysis.csv"))
matplot(realdata$t, realdata[,-1], type="b")

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = nrow(realdata),
    noise = c(NaN, 0.005, NaN, 0.005),
    kernel = "generalMatern",
    seed = 123,
    bandsize = 100,
    hmcSteps = 100,
    n.iter = 20001,
    stepSizeFactor = 0.01,
    linfillspace = 0.5, 
    t.end = 70,
    modelName = "Michaelis-Menten"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(5, 12.0341687, 2.03506154),
  x0 = c(0.1, 1, 0, 0),
  phi = cbind(c(0.2, 50), c(1, 50), c(1, 50), c(1, 50)),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(3,5)], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = round(realdata$t / config$linfillspace) * config$linfillspace)
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
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
  fOde=magi:::MichaelisMentenModelODE,
  fOdeDx=magi:::MichaelisMentenModelDx,
  fOdeDtheta=magi:::MichaelisMentenModelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="Michaelis-Menten"
)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- pram.true$sigma
sigma_fixed[1] <- 0.005
sigma_fixed[3] <-  0.005
sigma_fixed[2] <- 1e-4
sigma_fixed[4] <- 1e-4

xsim[1,-1] <- pram.true$x0

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control = 
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = config$stepSizeFactor,
                                  burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, discardBurnin=TRUE, useFixedSigma=TRUE))

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
pram.true$sigma[1] <- pram.true$sigma[3] <- 0
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
tail(gpode$theta)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)

#' TODO on simulated data
#' repeated experiments

# real-data ----
xreal <- xsim
xreal[1,-1] <- NA
colnames(xreal) <- c("time", "E", "S", "ES", "P")
xreal$S[is.finite(xreal$S)] <- realdata$S
xreal$P[is.finite(xreal$P)] <- realdata$P

xreal$ES[1] <- 0
xreal$P[1] <- 0
xreal$E[1] <- 0.1
xreal$S[1] <- 1

sigma_fixed[1] <- 1e-4  # cannot be too small due to numerical stability coded in c++ 
sigma_fixed[3] <-  1e-4  # 1e-4 is the smallest value

pram.true$phi = cbind(c(1, 50), c(1, 50), c(1, 50), c(0.2, 50))

result <- magi::MagiSolver(xreal[,-1], dynamicalModelList, xreal$time, control = 
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = config$stepSizeFactor,
                                  burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, discardBurnin=TRUE, useFixedSigma=TRUE))

gpode <- result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, dynamicalModelList$fOde(theta[t,], xsampled[t,,], xreal$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = dynamicalModelList$fOde(pram.true$theta, data.matrix(xtrue[,-1]), xtrue$time)

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xreal)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}

gpode$lglik <- gpode$lp
pram.true$sigma[1] <- pram.true$sigma[3] <- 0
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "-real-data.pdf"),
  xtrue, dotxtrue, xreal, gpode, pram.true, config, odemodel)
tail(gpode$theta)

tail(gpode$sigma)
apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)


#' TODO on real-data
#' MLE benchmark
#' still not sure about the "correct" initial condition, try checking with other papers.
