#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 15,
    noise = rep(0.001, 5), # 0.001 = low noise, 0.01 = high noise
    kernel = "generalMatern",
    seed = 1365546660, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    #loglikflag = "withmeanBand",
    #bandsize = 40,
    #hmcSteps = 100,
    n.iter = 300001,
    burninRatio = 0.50,
    #stepSizeFactor = 0.06,
    #filllevel = 3,
    #linfillspace = 0.1,
    t.end = 100,
    modelName = "PTrans-Dondel"
    #temperPrior = TRUE,
    #useFrequencyBasedPrior = TRUE,
    #useScalerSigma = FALSE,
    #useFixedSigma = FALSE,
    #linearizexInit = TRUE,
    #useExoSigma = TRUE,
    #max.epoch = 1
  )
}

# initialize global parameters, true x, simulated x ----------------------------
pram.true <- list(
  theta=c(0.07, 0.6,0.05,0.3,0.017,0.3),
  x0 = c(1,0,1,0,0),
  sigma=config$noise
)

times <- seq(0,100,length=1001)

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

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

ptransmodel <- list(
  fOde=gpds:::ptransmodelODE,
  fOdeDx=gpds:::ptransmodelDx,
  fOdeDtheta=gpds:::ptransmodelDtheta,
  thetaLowerBound=rep(0,6),
  thetaUpperBound=rep(100,6)
)


## Dondelinger demo on PTrans.

library(deSolve)
library(deGradInfer)

VG_func <- function(t, X, params) {
  MM <- ((params[5]*X[,5])/(params[6]+X[,5]))# Michaelis-Menten term
  dxdt <- cbind( - params[1]*X[,1] - params[2]*X[,1]*X[,3] +params[3]*X[,4],# S
                params[1]*X[,1],# dS
                - params[2]*X[,1]*X[,3] + params[3]*X[,4] +MM,# R
                params[2]*X[,1]*X[,3] - params[3]*X[,4] -params[4]*X[,4],# RS
                params[4]*X[,4] - MM# Rpp
        )
        return(dxdt)
}

# Generate data
timeTest <- xsim[,1]  #c(0,1,2,4,5,7,10,15,20,30,40,50,60,80,100)
dataTest <- as.matrix(xsim[,2:ncol(xsim)])

## noise 0.01 case
#dataTest = ode(c(1,0,1,0,0), timeTest,function(t,y,params)list(VG_func(t,matrix(y,1,length(y)),params)),c(0.07,0.6,0.05,0.3,0.017,0.3))
#dataTest = dataTest[,2:6] +rnorm(dim(dataTest)[1]*5,0,0.01)
agm.result =agm(data=dataTest,time=timeTest,ode.system=VG_func, numberOfParameters=length(pram.true$theta), noise.sd = config$noise[1], maxIterations = config$n.iter, showProgress = TRUE)
### you have to input noise.sd, method does not sample it


#### Use our plotting codes
config$n.iter <- config$n.iter / 25

xsampled <- (sapply(agm.result$x.samples, function(x) x, simplify="array"))
thetasampled <- agm.result$posterior.samples
lglik <- agm.result$ll

burnin <- as.integer(config$n.iter*config$burninRatio)
gpode <- list(theta= thetasampled[-(1:burnin),],
              xsampled= xsampled[-(1:burnin),,],
              lglik=  lglik[-(1:burnin)],
              sigma=  matrix(0, nrow = config$n.iter-burnin, ncol = ncol(xsim)-1)) # not sampled in this method
gpode$fode <- sapply(1:length(gpode$lglik), function(t) 
  with(gpode, gpds:::ptransmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
gpode$fode <- aperm(gpode$fode, c(3,1,2))

dotxtrue = gpds:::ptransmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

odemodel <- list(times=times, modelODE=modelODE, xtrue=xtrue)

outDir <- "./"

gpds:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".pdf"), 
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)
