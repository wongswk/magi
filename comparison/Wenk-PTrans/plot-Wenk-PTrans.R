# Use our plotting codes for Wenk output
library(gpds)
dataDir <- '/home/s246wong/github/FGPGM-output/noise0.001/'

config <- list(
  nobs = 15,
  noise = rep(0.01, 5), # 0.001 = low noise, 0.01 = high noise
  seed = 12345,   #### not used just placeholder
  n.iter = 300001,  ## change to iterations used in Wenk sampler
  burninRatio = 0.5,
  t.end = 100,
  modelName = "PTrans-Wenk"
)

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::ptransmodelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

xsim.obs <- cbind(xsim[,"time"], read.table(paste0(dataDir,"observations.csv")))   #### read Wenk observations
colnames(xsim.obs) <- c("time", 1:5)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

ptransmodel <- list(
  fOde=gpds:::ptransmodelODE,
  fOdeDx=gpds:::ptransmodelDx,
  fOdeDtheta=gpds:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(100,6)
)

samplesCpp <- as.matrix(read.table(paste0(dataDir, "MCMCMatrix.csv")))
xMean <- as.vector(as.matrix(read.table(paste0(dataDir, "meanMatrix.csv") )))
xSD <- as.vector(as.matrix(read.table(paste0(dataDir, "stdMatrix.csv") )))
thetaMag <- as.vector(as.matrix(read.table(paste0(dataDir, "thetaMagnitudes.csv") )))

llikId <- 0  ### llik not currently saved
xId <- (max(llikId)+1):(max(llikId)+length(data.matrix(xsim[,-1])))
thetaId <- (max(xId)+1):(max(xId)+length(ptransmodel$thetaLowerBound))
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
  with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- dataDir

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)


