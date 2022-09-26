library(magi)

outDir <- "../results/lac-operon/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) > 0){
  seed <- args
}else{
  seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9
}

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = c(0.01, 0.1, 0.1, 0.1, 0.001, 0.01, 0.01, 0.01, 0.01, 0.1),
    kernel = "generalMatern",
    seed = seed,
    bandsize = 100,
    hmcSteps = 100,
    niterHmc = 20001,
    stepSizeFactor = 0.001,
    filllevel = 0,
    t.end = 1200,
    modelName = "lac-operon"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(1, 0.02, 0.1, 0.005, 0.1, 1, 0.01, 0.1, 0.01, 0.03, 0.1, 0.001, 0.01, 0.002, 0.002, 0.01, 0.001),
  x0 = c(0, 50, 1000, 0, 1, 0, 100, 0, 0, 0),
  phi = cbind(c(1, 50), c(1, 50), c(1, 50), c(0.2, 50)),
  sigma=config$noise
)


times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::lacOperonODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(-1, -4)], type="l", lty=1)

# Compare to figure 3(a) in Barbuti, R., Gori, R., Milazzo, P., and Nasti, L. (2020). A survey of gene regula- tory networks modelling methods: from differential equations, to boolean and qualitative bioinspired models. Journal of Membrane Computing, 2(3):207– 226.
plot(xtrue[, "time"], xtrue[, "X10"], type="l", lty=1, main="Z")

matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,config$t.end,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim.obs <- xsim.obs[-1,] # remove first observation b.c. the rapid changes
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- setDiscretization(xsim.obs,config$filllevel)

dynamicalModelList <- list(
  fOde=magi:::lacOperonODE,
  fOdeDx=magi:::lacOperonDx,
  fOdeDtheta=magi:::lacOperonDtheta,
  thetaLowerBound=rep(0, 17),
  thetaUpperBound=rep(Inf, 17),
  name="lac-operon"
)

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList", xInitExogenous, pram.true$theta, xsim$time)


#' manually override estimated hyper-parameters for some components
#' GP smoothing gives bad result for rapidly decreasing curve
phiExogenous <- cbind(
  c(2.5, 600),
  c(100, 140),
  c(1000, 200),
  c(60, 100),
  c(0.01, 800),  # remove first obs
  c(1, 200),
  c(1, 300),
  c(0.5, 300),
  c(1, 400),
  c(35, 1000)
)
sigmaInit <- config$noise


#' works well except for theta0 which is the [i] component
#' even assume component [i] is known constant, k5 and k7 are still biased, other parameter inferences are good
#' trajectory RMSE are improved with known component [i]
#'
OursStartTime <- proc.time()[3] 
result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(xInit=xInitExogenous, niterHmc=config$niterHmc, stepSizeFactor = config$stepSizeFactor, phi=phiExogenous, sigma=sigmaInit, useFixedSigma=TRUE))
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
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-totalnoise",sum(config$noise), "xinitlin.pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

# TODO repeated experiments to get summary table, including coverage, trajectory RMSE, parameter RMSE
# tune the noise level and make the noise in different components comparable in scale

save.image(paste0(outDir, config$modelName,"-",config$seed, "-totalnoise",sum(config$noise),".rda"))

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-totalnoise",sum(config$noise),"-inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-totalnoise",sum(config$noise),"-inferred_theta.csv"))
