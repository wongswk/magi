library(magi)


outDir <- "../results/Michaelis-Menten-Va/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0("../results/Michaelis-Menten/", "hydrolysis.csv"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0){
  args <- as.numeric(args)
  seed <- args %% 100 + 1
  args <- args %/% 100
  
  noise_case = args %% 5
  args <- args %/% 5
  if(noise_case == 0){
    noise_scalar = 0.02 
  }else if(noise_case == 1){
    noise_scalar = 0.01
  }else if(noise_case == 2){
    noise_scalar = 0.005
  }else if(noise_case == 3){
    noise_scalar = 0.002
  }else if(noise_case == 4){
    noise_scalar = 0.001
  }
  noise = c(NA, noise_scalar, noise_scalar)
  
  # data_source_case <- args %% 2
  # args <- args %/% 2
  # if(data_source_case == 0){
  #   obs_source = "va-csv"
  # }else{
  #   obs_source = "vb-csv"
  # }
  obs_source = "vb-csv"
  
  hold_out_size <- args %% 15
  
  phi2 = 70
  phi = cbind(c(0.1, phi2), c(1, 30), c(1, 30))
  
  linfillspace = c(0.5)
  linfillcut = NULL
  phi_change_time = 0
  time_acce_factor = 1
  obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
  
  obs_time = c(0.5, 1, 2.5, 3.5, 4.5, 5.5, 7, 8.5, 9.5, 11, 12, 13.5, 15, 
               16, 18, 20, 21.5, 24, 27, 29.5, 32.5, 35.5, 39.5, 45, 55, 69)
  
  obs_time = obs_time[obs_keep]
  t.truncate = obs_time[length(obs_time) - hold_out_size]

}else{
  noise_scalar <- 0.01
  seed <- 1
  noise = c(NA, noise_scalar, noise_scalar)
  linfillspace = c(0.5)
  linfillcut = NULL
  phi = cbind(c(0.1, 70), c(1, 30), c(1, 30))
  phi_change_time = 0
  time_acce_factor = 1
  obs_keep = 1:26
  obs_source = "vb-csv"
  t.truncate = 70
}


# set up configuration if not already exist ------------------------------------

config <- list(
  nobs = nrow(realdata),
  noise = noise,
  kernel = "generalMatern",
  seed = seed,
  bandsize = 40,
  hmcSteps = 100,
  n.iter = 8001,
  linfillspace = linfillspace, 
  linfillcut = linfillcut,
  t.end = 70,
  t.start = 0,
  obs_start_time = 0,
  phi_change_time = phi_change_time,
  time_acce_factor = time_acce_factor,
  t.truncate = t.truncate,
  obs_keep = obs_keep,
  useMean = TRUE,
  phi = phi,
  skip_visualization = TRUE,
  obs_source = obs_source,
  modelName = "Michaelis-Menten-Va"
)

if(is.null(config$skip_visualization)){
  matplot(realdata$t, realdata[,-1], type="b")
}

# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list(
  theta=c(0.35, 0.2, 2.54),
  # x0 = c(0.08277011, 0.7500533, 0.2327168),
  x0 = c(0.1, 1, 0),
  phi = config$phi,
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelVaODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
if(is.null(config$skip_visualization)){
  matplot(xtrue[, "time"], xtrue[, c(3,4)], type="l", lty=1)  
  matplot(realdata$t, realdata[,-1], type="p", add=TRUE)
  matplot(realdata$t, realdata[,-1]/2, type="p", add=TRUE)
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
colnames(xsim.obs)[-1] <- c("E", "S", "P")
xsim.obs$E <- NULL
# write.csv(xsim.obs, paste0(outDir, "/va_xsim_obs_seed", config$seed, ".csv"))

if(config$obs_source == "vb-csv"){
  xsim.obs <- read.csv(paste0("../results/Michaelis-Menten-Vb4p/noise",noise_scalar,"/vb_xsim_obs_seed",config$seed,".csv"), row.names=1)
  xsim.obs$E <- NaN
  xsim.obs <- xsim.obs[,c("time", "E", "S", "P")]
  print("obs_source 'vb-csv'")
}else if(config$obs_source == "va-csv"){
  xsim.obs <- read.csv(paste0("../results/Michaelis-Menten-Va/noise",noise_scalar,"/va_xsim_obs_seed",config$seed,".csv"), row.names=1)
  xsim.obs$E <- NaN
  xsim.obs <- xsim.obs[,c("time", "E", "S", "P")]
  print("obs_source 'va-csv'")
}else if(config$obs_source == "va-sim"){
  xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
  colnames(xsim.obs)[-1] <- c("E", "S", "P")
  print("obs_source 'va-sim'")
}else{
  stop("obs_source not correct")
}

xsim.obs <- xsim.obs[config$obs_keep,]
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
xsim[xsim$time > config$t.truncate, -1] <- NaN

xsim.obs <- xsim.obs[xsim.obs$time >= config$t.start, ]
xsim <- xsim[xsim$time >= config$t.start, ]


# cpp inference ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenModelVaODE,
  fOdeDx=magi:::MichaelisMentenModelVaDx,
  fOdeDtheta=magi:::MichaelisMentenModelVaDtheta,
  thetaLowerBound=c(0,-100,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="Michaelis-Menten-Va"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

config$ndis <- nrow(xsim)

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

# MAGI off-the-shelf ----
# sampler with a good phi supplied, no missing component
# hyper-parameters affect the inference of missing components (especially if initial condition is not known

# xInitExogenous <- data.matrix(xsim[,-1])
# for (j in c(2,3)){
#   xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
#   idx <- which(is.na(xInitExogenous[, j]))
#   xInitExogenous[idx, j] <- xInitExogenous[idx[1] - 1, j]
# }
# xInitExogenous[-1, 1] <- 0.1

# xInitExogenous <- sapply(xtrueFunc, function(f) f(xsim$time))
xInitExogenous <- NULL

stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
# if(config$t.start == 0){
for(j in 1:3){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0
  }
}
# }


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
config$phi <- NULL

gpode$lglik <- gpode$lp
pram.true$sigma <- sigma_fixed
gpode$theta <- cbind(gpode$theta, (gpode$theta[,2]+gpode$theta[,3])/gpode$theta[,1])
pram.true$theta <- c(pram.true$theta, (pram.true$theta[2]+pram.true$theta[3])/pram.true$theta[1])

if(!is.null(config$linfillcut)){
  config$linfillcut <- paste(round(config$linfillcut, 2), collapse = ";")
  config$linfillspace <- paste(round(config$linfillspace, 2), collapse = ";")
}else{
  config$linfillcut <- NULL
}
config$obs_keep <- paste(c(config$obs_keep[1:5], ".."), collapse = ";")

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-fill", config$linfillspace,"-noise",
         sum(config$noise, na.rm = TRUE), "-phi", sum(pram.true$phi),
         "-data", config$obs_source,
         "-time", config$t.start,"to", config$t.truncate,
         "-obs_keep", config$obs_keep,
         "-linfillcut", config$linfillcut,
         "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
         ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
tail(gpode$theta)

apply(gpode$xsampled[,1,], 2, median)
apply(gpode$theta, 2, median)

save.image(paste0(outDir, config$modelName,"-",config$seed,"-fill", config$linfillspace,"-noise",
                  sum(config$noise, na.rm = TRUE), "-phi", sum(pram.true$phi),
                  "-data", config$obs_source,
                  "-time", config$t.start,"to", config$t.truncate,
                  "-obs_keep", config$obs_keep,
                  "-linfillcut", config$linfillcut,
                  "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
                  ".rda"))
