library(magi)

outDir <- "../results/lac-operon/"

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = rep(0.001, 10),
    kernel = "generalMatern",
    seed = 123,
    bandsize = 100,
    hmcSteps = 100,
    niterHmc = 2001,
    stepSizeFactor = 0.01,
    filllevel = 1,
    t.end = 1000,
    modelName = "lac-operon"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(1, 0.02, 0.1, 0.005, 0.1, 1, 0.01, 0.1, 0.01, 0.03, 0.1, 0.01, 0.01, 0.002, 0.002, 0.01, 0.001),
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
matplot(xtrue[, "time"], xtrue[, -4], type="l", lty=1)

# FIXME different from figure 3(a) in Barbuti, R., Gori, R., Milazzo, P., and Nasti, L. (2020). A survey of gene regula- tory networks modelling methods: from differential equations, to boolean and qualitative bioinspired models. Journal of Membrane Computing, 2(3):207– 226.
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

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- setDiscretization(xsim.obs,config$filllevel)

dynamicalModelList <- list(
  fOde=magi:::lacOperonODE,
  fOdeDx=magi:::lacOperonDx,
  fOdeDtheta=magi:::lacOperonDtheta,
  thetaLowerBound=rep(0, 17),
  thetaUpperBound=rep(Inf, 17)
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
                            "generalMatern")
  phiExogenous[,j] <- hyperparam$phi
  sigmaInit[j] <- hyperparam$sigma
  plot(xsim.obs$time, xsim.obs[,j+1], main=paste0("component ", j))
  lines(xtrue$time, xtrue[,j+1], col=2)
  mtext(paste0("sigma = ", round(sigmaInit[j], 3), 
               "; phi = ", paste0(round(phiExogenous[,j], 3), collapse = ",")))
}


#' manually override estimated hyper-parameters for some components
#' GP smoothing gives bad result for rapidly decreasing curve
phiExogenous[1,10] <- 0.01


OursStartTime <- proc.time()[3] 
result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(xInit = xInitExogenous, niterHmc=config$niterHmc, stepSizeFactor = config$stepSizeFactor, phi=phiExogenous, sigma=sigmaInit))
OursTimeUsed <- proc.time()[3] - OursStartTime
