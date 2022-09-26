# Summarize the results
library(magi)

# remove results that don't have common seed
rdaDir <- "../results/repressilator-gene-regulation-log/"
subdirs <- c(
  "../results/repressilator-gene-regulation-log/"
)

# get the csv quick summary first ----
config <- list()
config$modelName <- "repressilator-gene-regulation-log"
config$filllevel <- 0

theta_csv <- paste0("partobs-fill",config$filllevel,"-inferred_theta.csv")
trajectory_csv <- paste0("partobs-fill",config$filllevel,"-inferred_trajectory.csv")

for (rdaDir in subdirs){
  all_files <- list.files(rdaDir)
  theta_csv_files <- all_files[grep(theta_csv, all_files)]
  trajectory_csv_files <- all_files[grep(trajectory_csv, all_files)]
  
  common_seed <- intersect(
    gsub(".*log-([0-9]+)-partobs.*", "\\1", theta_csv_files),
    gsub(".*log-([0-9]+)-partobs.*", "\\1", trajectory_csv_files)
  )
  setdiff(1:1000, as.numeric(common_seed)) 
  # 217 575 613 641 642 764
  
  inferred_theta_all <- sapply(common_seed, function(seed){
    read.csv(paste0(rdaDir, config$modelName,"-",seed,"-partobs-fill",config$filllevel,"-inferred_theta.csv"))$x
  })
  
  inferred_trajectory_all <- sapply(common_seed, function(seed){
    data.matrix(read.csv(paste0(rdaDir, config$modelName,"-",seed,"-partobs-fill",config$filllevel,"-inferred_trajectory.csv")))
  }, simplify = "array")
  
  print(apply(inferred_theta_all, 1, mean))
  print(apply(inferred_theta_all, 1, sd))
    
}

  
for (rdaDir in subdirs){
  rm(list=setdiff(ls(), c("rdaDir", "subdirs", "config")))
  rdaDir <- paste0(rdaDir, "/")
  rdaDirSummary <- rdaDir
  print(rdaDir)
  
  pdf_files <- list.files(rdaDir)
  rda_files <- pdf_files[grep(paste0(".*log-([0-9]+)-partobs-fill",config$filllevel,".*\\.rda"), pdf_files)]
  
  
  ## Helper function adapted from Visualization to extract trajectories and RMSE
  rmsePostSamples <- function(xtrue, dotxtrue, xsim, gpode, param, config, odemodel=NULL){
    xpostmean <- apply(gpode$xsampled, 2:3, mean)
    
    times <- sort(unique(round(c(odemodel$times, xsim$time, xtrue[,"time"]), 7)))
    xdesolveTRUE <- deSolve::ode(y = param$x0, times = times, func = odemodel$modelODE, parms = param$theta)
    
    mapId <- which.max(gpode$lglik)
    ttheta <- gpode$theta[mapId,]
    tx0 <- gpode$xsampled[mapId,1,]
    xdesolveMAP <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    ttheta <- colMeans(gpode$theta)
    tx0 <- colMeans(gpode$xsampled[,1,])
    xdesolvePM <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    rowId <- sapply(xsim$time, function(x) which(abs(x-times) < 1e-6))
    xdesolveTRUE.obs <- xdesolveTRUE[rowId,-1]
    xdesolveMAP.obs <- xdesolveMAP[rowId,-1]
    xdesolvePM.obs <- xdesolvePM[rowId,-1]
    
    rmseInferredTrajectory <- sqrt(apply((xpostmean - xdesolveTRUE.obs)^2, 2, mean, na.rm=TRUE))
    
    return( list(xdesolveTRUE = xdesolveTRUE, xdesolvePM = xdesolvePM, xpostmean = xpostmean, rmseInferredTrajectory=rmseInferredTrajectory))
    
  }
  
  #### Grab the data for each seed and save in a list
  ours <- list()
  oursPostX <- list()
  oursPostTheta <- list()
  oursPostExpX <- list()
  oursExpXdesolvePM <- list()
  
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
      
      xsampledexp <- exp(gpode$xsampled)
      oursPostExpX_f <- cbind(
        apply(xsampledexp, 2:3, mean),
        apply(xsampledexp, 2:3, function(x) quantile(x, 0.025)),
        apply(xsampledexp, 2:3, function(x) quantile(x, 0.975)),
        apply(xsampledexp, 2:3, median)
      )
      
      ttheta <- colMeans(gpode$theta)
      exptx0 <- colMeans(xsampledexp[,1,])
      xdesolvePM <- deSolve::ode(y = log(exptx0), times = times, func = odemodel$modelODE, parms = ttheta)
      oursExpXdesolvePM_f <- exp(xdesolvePM)
      list(
        ours_f=ours_f,
        oursPostX_f=oursPostX_f,
        oursPostTheta_f=oursPostTheta_f,
        oursPostExpX_f=oursPostExpX_f,
        oursExpXdesolvePM_f=oursExpXdesolvePM_f
      )
    }, error = function(e) {
      return(as.character(e))
    })
  }, mc.cores = 64)
  
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
    oursPostExpX[[f]] <- outStorage[[f]]$oursPostExpX_f
    oursExpXdesolvePM[[f]] <- outStorage[[f]]$oursExpXdesolvePM_f
  }
  
  
  print(paste0(rdaDir,"summary-partobs-fill",config$filllevel, ".rda"))
  save.image(file=paste0(rdaDir,"summary-partobs-fill",config$filllevel, ".rda"))
}
