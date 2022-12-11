library(magi)


outDir <- "../results/Michaelis-Menten-Va/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0("../results/Michaelis-Menten/", "hydrolysis.csv"))

args <- commandArgs(trailingOnly = TRUE)


for (noise_scalar in c(0.02, 0.01, 0.005, 0.002, 0.001)){
  outDir <- paste0("../results/Michaelis-Menten-Va/", "noise", noise_scalar, "/")
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  
  for (seed in 1:100){
    noise = c(NA, noise_scalar, noise_scalar)
    linfillspace = c(0.5)
    linfillcut = NULL
    phi = cbind(c(0.1, 70), c(1, 30), c(1, 30))
    phi_change_time = 0
    time_acce_factor = 1
    obs_keep = 1:26
    obs_source = "va-csv"
    t.truncate = 70
    
    
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
      theta=c(1.17, 0, 0.32),
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
    write.csv(xsim.obs, paste0(outDir, "/va_xsim_obs_seed", config$seed, ".csv"))
    
  }
}

old_csv <- read.csv("../results/MM-model-comparison-wrong-ODE/Michaelis-Menten-Va.csv")
xdesolveTRUE <- deSolve::ode(y = pram.true$x0, times = old_csv$time, func = modelODE, parms = pram.true$theta)
colnames(xdesolveTRUE) <- colnames(old_csv)[-1]
write.csv(xdesolveTRUE, "../results/Michaelis-Menten-Va.csv")
