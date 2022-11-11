library(magi)

outDir <- "../results/Michaelis-Menten-Vb4p/"
realdata <- read.csv(paste0("../results/Michaelis-Menten/", "hydrolysis.csv"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

if(!exists("externalFlag")){
  externalFlag <- FALSE
}

# set up configuration if not already exist ------------------------------------
if(!externalFlag){
  config <- list(
    nobs = nrow(realdata),
    noise = c(NaN, 0.01, NaN, 0.01),
    kernel = "generalMatern",
    seed = 123,
    bandsize = 100,
    hmcSteps = 100,
    n.iter = 5001,
    t.start = 0,
    obs_keep=setdiff(1:26, c(1,2,4,6,8,11)),
    linfillspace = c(0.2, 0.5), 
    linfillcut = 3,
    stepSizeFactor = 0.01,
    linfillspace = 0.5, 
    phi_change_time = 3,
    time_acce_factor = 4,
    t.end = 70,
    modelName = "Michaelis-Menten-Vb4p"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(0.636, 0.0, 10.8, 0.0),
  x0 = c(0.1, 1, 0, 0),
  phi = cbind(c(0.1, 180), c(1, 30), c(0.1, 70), c(0.5, 30))
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelVb4pODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
sum(log(xtrue[-1,]))
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
matplot(xtrue[, "time"], xtrue[, c(3,5)], type="l", lty=1)
matplot(realdata$t, realdata[,-1], type="p", add=TRUE)
matplot(realdata$t, realdata[,-1]/2, type="p", add=TRUE)


# theta <- pram.true$theta
# k1 = theta[1]
# kn1 = theta[2]
# k2 = theta[3]
# (k1 * k2) / (kn1 + k2)

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
xsim.obs <- xsim.obs[config$obs_keep,]
xsim.obs <- rbind(c(0, pram.true$x0), xsim.obs)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

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

xsim.obs <- xsim.obs[xsim.obs$time >= config$t.start, ]
xsim <- xsim[xsim$time >= config$t.start, ]

# cpp inference under model Vb4p ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenModelVb4pODE,
  fOdeDx=magi:::MichaelisMentenModelVb4pDx,
  fOdeDtheta=magi:::MichaelisMentenModelVb4pDtheta,
  thetaLowerBound=c(0,-100,0,-100),
  thetaUpperBound=c(Inf,Inf,Inf,Inf),
  name="Michaelis-Menten-Vb4p"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

config$ndis <- config$t.end / config$linfillspace + 1

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
for(j in 1:4){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0  
  }
}

# xInitExogenous <- data.matrix(xsim[,-1])
# for (j in c(2,3)){
#   xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
#   idx <- which(is.na(xInitExogenous[, j]))
#   xInitExogenous[idx, j] <- xInitExogenous[idx[1] - 1, j]
# }
# xInitExogenous[-1, 1] <- 0.1

xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
# xInitExogenous <- NULL

distSignedCube <- array(NA, dim=c(nrow(xsim), nrow(xsim), ncol(xsim)-1))
for(j in 1:(ncol(xsim)-1)){
  for(i in 1:nrow(xsim)){
    distSignedCube[,i,j] = xsim$time - xsim$time[i]  
  }
}
tvec_accelarated = xsim$time
tvec_accelarated = tvec_accelarated - config$phi_change_time
tvec_accelarated[tvec_accelarated < 0] = tvec_accelarated[tvec_accelarated < 0] * config$time_acce_factor
for(j in c(1)){
  for(i in 1:nrow(xsim)){
    distSignedCube[,i,j] = tvec_accelarated - tvec_accelarated[i]  
  }
}

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control =
                             list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = stepSizeFactor,
                                  xInit = xInitExogenous, burninRatio = 0.5, phi = pram.true$phi, sigma=sigma_fixed, discardBurnin=TRUE, useFixedSigma=TRUE,
                                  skipMissingComponentOptimization=TRUE, useMean=config$useMean, useBand=FALSE, priorTemperature=NULL, distSignedCube=distSignedCube))


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
if(!is.null(config$linfillcut)){
  config$linfillcut <- paste(round(config$linfillcut, 2), collapse = ";")
  config$linfillspace <- paste(round(config$linfillspace, 2), collapse = ";")
}else{
  config$linfillcut <- NULL
}
config$obs_keep <- paste(c(config$obs_keep[1:5], ".."), collapse = ";")

gpode$lglik <- gpode$lp
pram.true$sigma <- sigma_fixed

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-fill", config$linfillspace,"-noise", 
         sum(config$noise, na.rm = TRUE), "-phi", sum(pram.true$phi),"-useMean", config$useMean,
         "-time", config$t.start,"to", config$t.truncate,"obsstart",config$obs_start_time, 
         "-obs_keep", config$obs_keep,
         "-linfillcut", config$linfillcut,
         "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
         ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)

#' #' TODO 
#' #' on simulated data repeated experiments
#' 
#' # use model A to fit model B data ----
#' 
#' dynamicalModelVa <- list(
#'   fOde=magi:::MichaelisMentenModelVaODE,
#'   fOdeDx=magi:::MichaelisMentenModelVaDx,
#'   fOdeDtheta=magi:::MichaelisMentenModelVaDtheta,
#'   thetaLowerBound=c(0,0,0),
#'   thetaUpperBound=c(Inf,Inf,Inf),
#'   name="Michaelis-Menten-Va"
#' )
#' 
#' xsim_va <- xsim[,c(1,2,3,4,6)]
#' 
#' sigma_va <- sigma_fixed
#' sigma_va <- sigma_va[c(1,2,3,5)]
#' 
#' phi_va <- pram.true$phi[,c(1,2,3,5)]
#' 
#' # MAGI off-the-shelf ----
#' # sampler with a good phi supplied, no missing component
#' # hyper-parameters affect the inference of missing components (especially if initial condition is not known
#' 
#' OursStartTime <- proc.time()[3]
#' 
#' result <- magi::MagiSolver(xsim_va[,-1], dynamicalModelVa, xsim_va$time, control = 
#'                              list(bandsize=config$bandsize, niterHmc=config$n.iter, nstepsHmc=config$hmcSteps, stepSizeFactor = config$stepSizeFactor,
#'                                   burninRatio = 0.5, phi = phi_va, sigma=sigma_va, discardBurnin=FALSE, useFixedSigma=TRUE,
#'                                   skipMissingComponentOptimization=TRUE))
#' 
#' OursTimeUsed <- proc.time()[3] - OursStartTime
#' 
#' 
#' gpode <- result
#' gpode$fode <- sapply(1:length(gpode$lp), function(t)
#'   with(gpode, dynamicalModelVa$fOde(theta[t,], xsampled[t,,], xsim_va$time)), simplify = "array")
#' gpode$fode <- aperm(gpode$fode, c(3,1,2))
#' 
#' 
#' xtrue_va <- xtrue[,c(1,2,3,4,6)]
#' 
#' dotxtrue_va = dotxtrue[,c(1,2,3,4,6)]
#' 
#' modelODEVa <- function(t, state, parameters) {
#'   list(as.vector(magi:::MichaelisMentenModelVaODE(parameters, t(state), t)))
#' }
#' 
#' odemodel <- list(times=times, modelODE=modelODEVa, xtrue=xtrue)
#' 
#' for(j in 1:(ncol(xsim)-1)){
#'   config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
#' }
#' 
#' gpode$lglik <- gpode$lp
#' pram.true$sigma <- sigma_fixed
#' 
#' magi:::plotPostSamplesFlex(
#'   paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"),
#'   xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
#' tail(gpode$theta)
#' 
#' apply(gpode$xsampled[,1,], 2, median)
#' apply(gpode$theta, 2, median)
