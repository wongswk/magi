library(magi)

outDir <- "../results/fn/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = c(0.2, 0.2),
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    niterHmc = 20000,
    filllevel = 2,
    t.end = 20,
    modelName = "FN"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::fnmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,20,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- setDiscretization(xsim.obs,config$filllevel)

fnmodel <- list(
  fOde=magi:::fnmodelODE,
  fOdeDx=magi:::fnmodelDx,
  fOdeDtheta=magi:::fnmodelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="FN"
)

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], fnmodel, xsim$time, control = list(xInit = xInitExogenous, niterHmc=config$niterHmc, stepSizeFactor = 0.06))

OursTimeUsed <- proc.time()[3] - OursStartTime

gpode <- result
gpode$fode <- sapply(1:length(gpode$lp), function(t) 
  with(gpode, magi:::fnmodelODE(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = magi:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}

gpode$lglik <- gpode$lp
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-FN_inferred_theta.csv"))
