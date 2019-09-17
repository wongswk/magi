#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 41,
    noise = c(0.2, 0.2),  # Dondel must have known scalar noise 
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 100000,
    burninRatio = 0.50,
    stepSizeFactor = 0.06,
    filllevel = 2,
    modelName = "FN",
    temperPrior = TRUE,
    useFrequencyBasedPrior = TRUE,
    useScalerSigma = FALSE,
    useFixedSigma = FALSE,
    npostplot = 10
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
if(config$temperPrior){
  config$priorTemperature <- config$ndis / config$nobs  
}else{
  config$priorTemperature <- 1
}

if(config$loglikflag == "withmeanBand"){
  config$useMean = TRUE
  config$useBand = TRUE
}else if(config$loglikflag == "band"){
  config$useMean = FALSE
  config$useBand = TRUE
}else if(config$loglikflag == "withmean"){
  config$useMean = TRUE
  config$useBand = FALSE
}else if(config$loglikflag == "usual"){
  config$useMean = FALSE
  config$useBand = FALSE
}

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,20,length=241)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
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

xsim <- insertNaN(xsim.obs,config$filllevel)


## Dondelinger demo on FN

library(deGradInfer)
dataTest = data.matrix(xsim.obs[,-1])
timeTest = xsim.obs$time

FN_func <- function(t, x, params) {
  theta <- params
  V <- x[, 1]
  R <- x[, 2]
  result <- array(0, c(nrow(x), ncol(x)))
  result[, 1] = theta[3] * (V - V^3/3 + R)
  result[, 2] = -1/theta[3] * (V - theta[1] + theta[2] * R)
  result
}

agm.result = agm(data=dataTest,
                 time=timeTest,
                 ode.system=FN_func, 
                 numberOfParameters=3, 
                 noise.sd = config$noise[1],
                 maxIterations = config$n.iter, 
                 showProgress = TRUE)
### you have to input noise.sd, method does not sample it, for FN system we don't specify noise


#### Use our plotting codes
xsampled <- (sapply(agm.result$x.samples, function(x) x, simplify="array"))
thetasampled <- agm.result$posterior.samples
lglik <- agm.result$ll

config$n.iter <- config$n.iter / 25

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta= thetasampled[-(1:burnin),],
              xsampled= xsampled[-(1:burnin),,],
              lglik=  lglik[-(1:burnin)],
              sigma=  matrix(0, nrow = config$n.iter-burnin, ncol = ncol(xsim)-1)) # not sampled in this method
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- "./"

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-Dondel-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim.obs, gpode, pram.true, config, odemodel)
