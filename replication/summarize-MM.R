# Summarize the results
library(magi)

# remove results that don't have common seed
rdaDir <- "../results/Michaelis-Menten/"
subdirs <- c(
  "../results/Michaelis-Menten/"
)

# get the csv quick summary first ----
for(phi2 in c(70, 90, 120)){
for(scenario in 0:3){
  
phi = cbind(c(0.1, phi2), c(1, 30), c(1, 30))

linfillspace = c(0.5)
linfillcut = NULL
phi_change_time = 0
time_acce_factor = 1
noise = c(NA, 0.02, 0.02)

if(scenario == 0){
  obs_keep = 1:26
}else if (scenario == 1){
  obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
}else if (scenario == 2){
  obs_keep = 4:26
}else if (scenario == 3){
  obs_keep = seq(2, 26, 2)
}

config <- list(
  nobs = 26,
  noise = noise,  # noise = c(0.01, 0.01, 0.01, 0.01), for fully observed case
  kernel = "generalMatern",
  bandsize = 40,
  hmcSteps = 100,
  n.iter = 5001,
  linfillspace = linfillspace, 
  linfillcut = linfillcut,
  t.end = 70,
  t.start = 0,
  obs_start_time = 0,
  phi_change_time = phi_change_time,
  time_acce_factor = time_acce_factor,
  t.truncate = 70,
  obs_keep = obs_keep,
  useMean = TRUE,
  phi = phi,
  skip_visualization = TRUE,
  modelName = "Michaelis-Menten-Reduced"
)

rm(list=setdiff(ls(), c("rdaDir", "subdirs", "config", "scenario", "phi2")))
rdaDir <- paste0(rdaDir, "/")
rdaDirSummary <- rdaDir
print(rdaDir)

if(!is.null(config$linfillcut)){
  config$linfillcut <- paste(round(config$linfillcut, 2), collapse = ";")
  config$linfillspace <- paste(round(config$linfillspace, 2), collapse = ";")
}else{
  config$linfillcut <- NULL
}

config$obs_keep <- paste(c(config$obs_keep[1:5], ".."), collapse = ";")


pdf_files <- list.files(rdaDir)
rda_files <- pdf_files[grep(paste0(config$modelName,"-.*-fill", config$linfillspace,"-noise", 
                                   sum(config$noise, na.rm = TRUE), "-phi", sum(config$phi),"-useMean", config$useMean,
                                   "-time.*", "-obs_keep", config$obs_keep, "-linfillcut", config$linfillcut,
                                   "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
                                   "\\.rda$"), pdf_files)]
print(rda_files[1])
print(length(rda_files))

## Helper function adapted from Visualization to extract trajectories and RMSE
rmsePostSamples <- function(xtrue, dotxtrue, xsim, gpode, param, config, odemodel=NULL){
  xpostmean <- apply(gpode$xsampled, 2:3, mean)
  
  times <- sort(unique(round(c(odemodel$times, xsim$time, xtrue[,"time"]), 7)))
  xdesolveTRUE <- deSolve::ode(y = param$x0, times = times, func = odemodel$modelODE, parms = param$theta)
  
  mapId <- which.max(gpode$lglik)
  ttheta <- gpode$theta[mapId,]
  tx0 <- gpode$xsampled[mapId,1,]
  starttime <- xsim$time[which(is.finite(xsim[,2]))[1]]
  
  ttheta <- colMeans(gpode$theta)
  tx0 <- colMeans(gpode$xsampled[,1,])
  xdesolvePM <- deSolve::ode(y = tx0, times = times[times >= starttime], func = odemodel$modelODE, parms = ttheta)
  
  return( list(xdesolveTRUE = xdesolveTRUE, xdesolvePM = xdesolvePM, xpostmean = xpostmean))
  
}

#### Grab the data for each seed and save in a list
ours <- list()
oursPostX <- list()
oursPostTheta <- list()

library(parallel)

outStorage <- mclapply(rda_files, function(f){
  tryCatch({
    load(paste0(rdaDirSummary, f), envir = .GlobalEnv)
    if (!exists("param_restricted")){
      param_restricted <- pram.true
    }
    
    ours_f <- rmsePostSamples(xtrue, dotxtrue, xsim, gpode, param_restricted, config, odemodel)
    oursPostX_f <- cbind(
      apply(gpode$xsampled, 2:3, mean),
      apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.025)),
      apply(gpode$xsampled, 2:3, function(x) quantile(x, 0.975)),
      apply(gpode$xsampled, 2:3, median)
    )
    oursPostTheta_f <- cbind(
      apply(gpode$theta, 2, mean),
      apply(gpode$theta, 2, function(x) quantile(x, 0.025)),
      apply(gpode$theta, 2, function(x) quantile(x, 0.975)),
      apply(gpode$theta, 2, median)
    )
    
    list(
      ours_f=ours_f,
      oursPostX_f=oursPostX_f,
      oursPostTheta_f=oursPostTheta_f
    )
  }, error = function(e) {
    return(as.character(e))
  })
}, mc.cores = 16)

valid_result_id <- sapply(1:length(rda_files), function(f) length(outStorage[[f]]) > 1)
error_msg <- outStorage[!valid_result_id]
error_file <- rda_files[!valid_result_id]
outStorage <- outStorage[valid_result_id]
rda_files <- rda_files[valid_result_id]

rda_files <- as.character(unlist(rda_files))
for (f in 1:length(rda_files)) {
  ours[[f]] <- outStorage[[f]]$ours_f
  oursPostX[[f]] <- outStorage[[f]]$oursPostX_f
  oursPostTheta[[f]] <- outStorage[[f]]$oursPostTheta_f
}

summary_filename = paste0(rdaDir,"summary-fill", config$linfillspace,"-noise", 
                          sum(config$noise, na.rm = TRUE), "-phi", sum(config$phi),"-useMean", config$useMean,
                          "-obs_keep", config$obs_keep, "-linfillcut", config$linfillcut,
                          "-time_changepoint", config$phi_change_time, "factor", config$time_acce_factor,
                          ".rda")

print(summary_filename)

save.image(file=summary_filename)
}
}
