library(magi)

outDir <- "../results/repressilator-gene-regulation-log/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) > 0){
  seed <- args
}else{
  seed <- 1
}

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = rep(0.15, 6),
    kernel = "generalMatern",
    seed = seed,
    niterHmc = 10001,
    stepSizeFactor = 0.001,
    filllevel = 0,
    t.end = 300,
    modelName = "repressilator-gene-regulation-log"
  )
}

csv_path = paste0(outDir, config$modelName,"-",config$seed,"-partobs-fill",config$filllevel,"-inferred_theta.csv")
if(file.exists(csv_path)){
  out = scan(csv_path, what="character")
  if(length(out) > 0){
    quit(save="no")  
  }
}


# initialize global parameters, true x, simulated x ----------------------------
alpha <- 240 # obtain from Fig 1b in Elowitz and Leibler (2000)
KM <- 40     # scale factor only, to convert protein number to match Fig 1c in paper
pram.true <- list(
  theta=c(0.001*alpha, alpha, 2, 1/5),  # alpha0/alpha = 0.001
  # initial condition cannot be the same, otherwise the system degenerates to two-components 
  # -- all the m and all the p will be the same
  x0 = log(c(0.4, 20, 40, 0.01, 0.01, 0.01)),
  # phi = cbind(c(1, 50), c(1, 50), c(1, 50)),
  sigma=config$noise
)


times <- seq(0,config$t.end,length=1001)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::repressilatorGeneRegulationLogODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
# Plot proteins only (times KM factor), compare to Fig 1c (left panel) in Elowitz and Leibler (2000)
matplot(xtrue[, "time"], xtrue[, -(1:4)] * KM, type="l", lty=1)
matplot(xtrue[, "time"], exp(xtrue[, -(1:4)]) * KM, type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,config$t.end,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}


xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xtrue[, "time"], xtrue[, -(1:4)], type="l", lty=1)
matplot(xsim.obs$time, xsim.obs[,-(1:4)], type="p", col=1:(ncol(xsim)-1), pch=20, add=TRUE)

xsim <- setDiscretization(xsim.obs,config$filllevel)

dynamicalModelList <- list(
  fOde=magi:::repressilatorGeneRegulationLogODE,
  fOdeDx=magi:::repressilatorGeneRegulationLogDx,
  fOdeDtheta=magi:::repressilatorGeneRegulationLogDtheta,
  thetaLowerBound=rep(0, 4),
  thetaUpperBound=rep(Inf, 4),
  name="repressilator-gene-regulation-log"
)

xInitExogenous <- data.matrix(xsim[,-1])
for (j in 1:(ncol(xsim)-1)){
  xInitExogenous[, j] <- approx(xsim.obs$time, xsim.obs[,j+1], xsim$time)$y
}

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList", xInitExogenous[-1,], pram.true$theta, xsim$time[-1])

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

if(config$filllevel == 0){
  phiExogenous <- rbind(rep(6, 6), rep(10, 6))
}else{
  phiExogenous <- rbind(rep(6, 6), rep(5, 6))  
}


# remove first obs
xsim <- xsim[-1,]

# protein levels missing
xsim[,5:7] <- NA
xsim.obs[,5:7] <- NA
sigmaInit <- config$noise

OursStartTime <- proc.time()[3] 

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(niterHmc=config$niterHmc, stepSizeFactor = config$stepSizeFactor, phi=phiExogenous, sigma=sigmaInit, useFixedSigma=TRUE))
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

save.image(paste0(outDir, config$modelName,"-",config$seed, "-partobs-fill",config$filllevel,".rda"))

magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "-partobs-fill",config$filllevel,".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

write.csv(apply(gpode$xsampled, 2:3, mean), paste0(outDir, config$modelName,"-",config$seed,"-partobs-fill",config$filllevel,"-inferred_trajectory.csv"))
write.csv(apply(gpode$theta, 2, mean), paste0(outDir, config$modelName,"-",config$seed,"-partobs-fill",config$filllevel,"-inferred_theta.csv"))
