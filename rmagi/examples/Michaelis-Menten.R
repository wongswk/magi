library(magi)

outDir <- "../results/Michaelis-Menten/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)


# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 11,
    noise = c(0.01, 0.01),
    seed = 123,
    niterHmc = 200,
    filllevel = 2,
    t.end = 70,
    modelName = "Michaelis-Menten"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(30, 2, 0.005),
  x0 = c(2, 1, 0, 0),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelODE(parameters, t(state))))
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

model <- list(
  fOde=magi:::MichaelisMentenModelODE,
  fOdeDx=magi:::MichaelisMentenModelDx,
  fOdeDtheta=magi:::MichaelisMentenModelDtheta,
  thetaLowerBound=c(0,0,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="MichaelisMenten"
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "Michaelis-Menten model", x, theta))

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

OursStartTime <- proc.time()[3]

result <- magi::MagiSolver(xsim[,-1], model, xsim$time, control = list(xInit = xInitExogenous, niterHmc=config$niterHmc, stepSizeFactor = 0.06))

OursTimeUsed <- proc.time()[3] - OursStartTime

gpode <- result
gpode$fode <- sapply(1:length(gpode$lp), function(t)
  with(gpode, model$fOde(theta[t,], xsampled[t,,], xsim$time)), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = model$fOde(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

for(j in 1:(ncol(xsim)-1)){
  config[[paste0("phiD", j)]] <- paste(round(gpode$phi[,j], 2), collapse = "; ")
}

gpode$lglik <- gpode$lp
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save(xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel, OursTimeUsed, file= paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-MichaelisMenten_inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-MichaelisMenten_inferred_theta.csv"))
