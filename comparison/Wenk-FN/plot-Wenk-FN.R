# Use our plotting codes for Wenk output
library(gpds)
dataDir <- './'

config <- list(
  nobs = 41,
  noise = c(0.2, 0.2),  # Dondel must have known scalar noise 
  seed = 12345,   #### not used just placeholder
  n.iter = 300001,  ## change to iterations used in Wenk sampler
  burninRatio = 0.5,
  t.end = 10,
  modelName = "FN"
)

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  sigma=config$noise
)


times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0, config$t.end, length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

xsim.obs <- cbind(xsim[,"time"], read.table(paste0(dataDir,"observations.csv")))   #### read Wenk observations
colnames(xsim.obs) <- c("time", 1:2)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))

llikId <- 0  ### llik not currently saved
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+3)
#sigmaId <- (max(thetaId)+1):(max(thetaId)+ncol(xsim[,-1]))  ## no sigma sampled

## Un-standardize their X's and untransform thetas
samplesCpp[,xId] <- t(apply(samplesCpp[,xId],1, function(x) xSD*x + xMean))
samplesCpp[,thetaId] <- t(apply(samplesCpp[,thetaId],1, function(x) x * 10^thetaMag))

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta= samplesCpp[-(1:burnin), thetaId],
              xsampled=array(samplesCpp[-(1:burnin), xId],
                             dim=c(config$n.iter-burnin, nrow(xsim), ncol(xsim)-1)),
              lglik=  rep(0, config$n.iter - burnin), ###samplesCpp[llikId,-(1:burnin)],
              sigma=  matrix(0, nrow = config$n.iter - burnin, ncol = ncol(xsim)-1)) # t(samplesCpp[sigmaId, -(1:burnin), drop=FALSE]))
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- dataDir

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-wenk-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)


