# Summarize the results
library(magi)

# remove results that don't have common seed
rdaDir <- "../results/Michaelis-Menten/"
subdirs <- c(
  "../results/Michaelis-Menten/"
)

# get the csv quick summary first ----
config <- list()
config$modelName <- "Michaelis-Menten"
config$linfillspace <- 0.5
config$noise = c(0.002, 0.02, 0.002, 0.02)
config$phi = cbind(c(0.1, 70), c(1, 30), c(0.1, 70), c(1, 30))

rm(list=setdiff(ls(), c("rdaDir", "subdirs", "config")))
rdaDir <- paste0(rdaDir, "/")
rdaDirSummary <- rdaDir
print(rdaDir)

pdf_files <- list.files(rdaDir)
rda_files <- pdf_files[grep(paste0(config$modelName,"-.*-fill", sum(config$linfillspace),"-noise", sum(config$noise, na.rm = TRUE), "-phi", sum(config$phi), "\\.rda"), pdf_files)]

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

print(paste0(rdaDir,"summary-fill",sum(config$linfillspace),"-noise", sum(config$noise, na.rm = TRUE), ".rda"))

save.image(file=paste0(rdaDir,"summary-fill",sum(config$linfillspace),"-noise", sum(config$noise, na.rm = TRUE), ".rda"))
