library(magi)

outDir <- "../results/Michaelis-Menten/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

#realdata <- read.csv(paste0(outDir, "hydrolysis.csv"))
#matplot(realdata$t, realdata[,-1], type="b")

obs.times <- c(2.5, 4.5, 7, 9.5, 11, 13.5, 15, 16, 18, 20, 21.5, 24, 27, 29.5, 32.5, 35.5, 39.5, 45, 55, 69)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = length(obs.times),
    noise = c(NA, 0.02, 0.02), 
    kernel = "generalMatern",
    seed = 123,
    n.iter = 5001,
    linfillspace = 0.5, 
    t.start = 0,
    t.end = 70,  
    phi = cbind(c(0.1, 70), c(1, 30), c(1, 30)),
    modelName = "Michaelis-Menten-Reduced"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list(
  theta=c(0.9, 0.75, 2.54),
  x0 = c(0.1, 1, 0),
  phi = config$phi,
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenReducedODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, c(3,4)], type="l", lty=1)  

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = round(obs.times / config$linfillspace) * config$linfillspace)
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))
xtest <- xsim

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
xsim.obs <- rbind(c(0, pram.true$x0), xsim.obs)

if(is.null(config$skip_visualization)){
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
  matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)
}

## Linearly interpolate using fixed interval widths
xsim <- setDiscretization(xsim.obs, by = config$linfillspace)

# cpp inference ----------------------------
dynamicalModelList <- list(
  fOde=magi:::MichaelisMentenReducedODE,
  fOdeDx=magi:::MichaelisMentenReducedDx,
  fOdeDtheta=magi:::MichaelisMentenReducedDtheta,
  thetaLowerBound=c(0,-100,0),
  thetaUpperBound=c(Inf,Inf,Inf),
  name="Michaelis-Menten-Reduced"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList",
                   data.matrix(xtest[,-1]), pram.true$theta, xtest$time)

sigma_fixed <- config$noise
sigma_fixed[is.na(sigma_fixed)] <- 1e-4

# MAGI off-the-shelf ----
xInitExogenous <- data.matrix(xsim[,-1])
for (j in c(2,3)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
  idx <- which(is.na(xInitExogenous[, j]))
  xInitExogenous[idx, j] <- xInitExogenous[idx[1] - 1, j]
}
xInitExogenous[-1, 1] <- 0.1

stepSizeFactor <- rep(0.01, nrow(xsim)*length(pram.true$x0) + length(dynamicalModelList$thetaLowerBound) + length(pram.true$x0))
for(j in 1:3){
  for(incre in 1:1){
    stepSizeFactor[(j-1)*nrow(xsim) + incre] <- 0  
  }
}

OursStartTime <- proc.time()[3]

gpode <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, control =
                             list(niterHmc=config$n.iter, stepSizeFactor = stepSizeFactor,
                                  xInit = xInitExogenous, phi = pram.true$phi, sigma=sigma_fixed,useFixedSigma=TRUE,
                                  skipMissingComponentOptimization=TRUE, verbose=TRUE))

OursTimeUsed <- proc.time()[3] - OursStartTime

# Add KM to parameters as a function of k1, k_{-1}, k2
gpode$theta <- cbind(gpode$theta, (gpode$theta[,2]+gpode$theta[,3])/gpode$theta[,1])
pram.true$theta <- c(pram.true$theta, (pram.true$theta[2]+pram.true$theta[3])/pram.true$theta[1])

par.table <- function(res) {
  par.est <- apply(cbind(res$theta[,-c(1:2)]), 2,
                   function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
  colnames(par.est) <- c("k_cat", "KM")
  rownames(par.est) <- c("Mean", "2.5%", "97.5%")
  signif(par.est, 3)
}

par.table(gpode)

