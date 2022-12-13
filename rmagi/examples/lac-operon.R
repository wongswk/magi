library(magi)

noisefac <- 0.05
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    noise = c(0.0236194228362474, 0.984379168607761, 0.040561105491602, 0.19224800630503,
              0.000104872420893985, 0.395628293932429, 99.3980508395468, 0.0263189487352539,
              0.00631925202022132, 0.00036869187211123)*noisefac, # min of each component * noisefac
    kernel = "generalMatern",
    seed = 118,
    #seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    hmcSteps = 500,
    niterHmc = 50001,
    fillinterval = 15,
    t.start = 1,
    t.end = 1201,
    modelName = "lac-operon",
    obs.times = c(seq(1,361,by=15), seq(391,601,by=30), 901, 1201)
  )
}

config$nobs = length(config$obs.times)
outDir <- paste0("../results/lac-operon/noise", noisefac, "-nobs", config$nobs, "-uneven/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(1, 0.02, 0.1, 0.005, 0.1, 1, 0.01, 0.1, 0.01, 0.03, 0.1, 0.001, 0.01, 0.002, 0.002, 0.01, 0.001),
  x0 = c(0, 50, 1000, 0, 1, 0, 100, 0, 0, 0),
  phi = cbind(c(1, 50), c(1, 50), c(1, 50), c(0.2, 50)),
  sigma=config$noise
)

times <- seq(0,config$t.end,by = 0.01)

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

xsim <- data.frame(time = config$obs.times)
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}

xsim.obs <- xsim
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- setDiscretization(xsim.obs, by = config$fillinterval)

dynamicalModelList <- list(
  fOde=magi:::lacOperonODE,
  fOdeDx=magi:::lacOperonDx,
  fOdeDtheta=magi:::lacOperonDtheta,
  thetaLowerBound=rep(0, 17),
  thetaUpperBound=rep(Inf, 17),
  name="lac-operon"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList", data.matrix(xsim.obs[,-1]), pram.true$theta, xsim$time)

phiExogenous <- cbind(
  c(2.5, 600),
  c(100, 140),
  c(1000, 200),
  c(60, 100),
  c(0.01, 800),
  c(1, 200),
  c(500, 300),
  c(0.5, 300),
  c(1, 400),
  c(35, 1000)
)
sigmaInit <- config$noise

OursStartTime <- proc.time()[3] 
result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(niterHmc=config$niterHmc, nstepsHmc = config$hmcSteps, phi=phiExogenous, sigma=sigmaInit, useFixedSigma=TRUE, verbose=TRUE))
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
  paste0(outDir, config$modelName,"-",config$seed,"-nobs",config$nobs,"-noise", config$noise[1], "xinitlin.pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)


par.table <- function(res) {
  par.est <- apply(cbind(res$theta[,-1]), 2,
                   function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))
  colnames(par.est) <- paste0('k', 1:16)
  rownames(par.est) <- c("Mean", "2.5%", "97.5%")
  signif(par.est, 3)
}

par.table(gpode)

tvecsolve <- seq(config$t.start,config$t.end,by = 0.1)
calcTraj <- function(res) {
  x0.est <- apply(res$xsampled[,1,],2,mean)
  theta.est <- apply(res$theta,2,mean)
  
  x <- deSolve::ode(y = x0.est, times = tvecsolve,
                    func = modelODE, parms = theta.est)
  x
}

recon <- calcTraj(gpode)
recon.obs <- subset(recon, time %in% xsim$time)
xtrue.obs <- subset(xtrue, time %in% xsim$time)[,-1]

# Trajectory RMSE
sqrt(colMeans((recon.obs - xtrue.obs)^2))

pdf(paste0(outDir, config$modelName,"-",config$seed,"-nobs",config$nobs,"-noise", config$noise[1], "-reconstruct.pdf"), width=12, height=6)
par(mfrow=c(2,5))
for (i in 1:10) {
  plot(c(min(tvecsolve), max(tvecsolve)), c(min(c(recon[,i+1], xtrue[xtrue$time >=1,i+1], xsim.obs[,i+1])), max(c(recon[,i+1], xtrue[xtrue$time>=1,i+1], xsim.obs[,i+1]))), type = 'n', ylab='', xlab=paste('Component', i))
  points(xsim[,1], xsim[,i+1], col='gray50')
  lines(xtrue$time[xtrue$time >=1], xtrue[xtrue$time >=1,i+1], col="red")
  lines(tvecsolve, recon[,i+1])
  
}
dev.off()

save.image(paste0(outDir, config$modelName,"-",config$seed,"-nobs",config$nobs,"-noise", config$noise[1], ".rda"))


# TODO repeated experiments to get summary table, including coverage, trajectory RMSE, parameter RMSE
# tune the noise level and make the noise in different components comparable in scale
