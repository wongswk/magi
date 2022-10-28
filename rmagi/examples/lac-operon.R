library(magi)

outDir <- "../results/lac-operon/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = c(0.0236194228362474, 0.984379168607761, 0.040561105491602, 0.19224800630503, 
      0.000104872420893985, 0.395628293932429, 99.3980508395468, 0.0263189487352539, 
      0.00631925202022132, 0.00036869187211123)*0.05,  # min of each component * 0.05
    kernel = "generalMatern",
    seed = 123,
    bandsize = 100,
    hmcSteps = 200,
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


times <- seq(0,config$t.end,length=100001)

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

phiExogenous <- matrix(0, nrow=2, ncol=ncol(xsim)-1)
sigmaInit <- rep(0, ncol(xsim)-1)
for (j in 1:(ncol(xsim)-1)){
  hyperparam <- gpsmoothing(xsim.obs[,j+1],
                            xsim.obs$time,
                            "generalMatern",
                            sigma = config$noise[j])
  phiExogenous[,j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
  plot(xsim.obs$time, xsim.obs[,j+1], main=paste0("component ", j))
  lines(xtrue$time, xtrue[,j+1], col=2)
  mtext(paste0("sigma = ", round(sigmaInit[j], 3),
               "; phi = ", paste0(round(phiExogenous[,j], 3), collapse = ", ")))
}


#' manually override estimated hyper-parameters for some components
#' GP smoothing gives bad result for rapidly decreasing curve
phiExogenous <- cbind(
  c(2.5, 600),
  c(100, 140),
  c(1000, 200),
  c(60, 100),
  c(0.01, 800),  # remove first obs
  c(1, 200),
  c(500, 300),  # increased phi1 to 500 
  c(0.5, 300),
  c(1, 400),
  c(35, 1000)
)
sigmaInit <- config$noise


#' assume component [i] is known constant, parameter inferences are good
#' trajectory RMSE are improved with known component [i]
#'
OursStartTime <- proc.time()[3] 
result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(niterHmc=config$niterHmc, nstepsHmc = config$hmcSteps, phi=phiExogenous, sigma=sigmaInit, useFixedSigma=TRUE))
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
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "xinitlin.pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

# TODO repeated experiments to get summary table, including coverage, trajectory RMSE, parameter RMSE
# tune the noise level and make the noise in different components comparable in scale
