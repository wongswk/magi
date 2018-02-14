rm(list=ls())
#### run with priorTempered phase 1 --------------------------------------------
library(gpds)
if(!exists("config")){
  config <- list(
    nobs = 11,
    noise = c(4,1,8)*0.2,
    kernel = "generalMatern",
    seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    loglikflag = "withmeanBand",
    bandsize = 20,
    hmcSteps = 500,
    n.iter = 5000,
    burninRatio = 0.50,
    stepSizeFactor = 1,
    filllevel = 3,
    modelName = "hes1"
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1
config$priorTemperature <- config$ndis / config$nobs
if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
}else{
  baseDir <- "~/Workspace/DynamicSys/results/"  
}
outDir <- with(config, paste0(baseDir, modelName, "-", loglikflag,"-", kernel,
                              "-nobs",nobs,"-noise", paste(round(noise,3), collapse = "_"),
                              "-ndis",ndis,"/"))
system(paste("mkdir -p", outDir))

pram.true <- list(
  theta = c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3),
  x0 = c(0.5, 2, 1),
  phi = c(122.4027613, 41.8511396,  
          56.5612956, 91.4189948,
          164.3556832, 11.9474091),
  sigma = config$noise
)
times <- seq(0, 60*4, by = 0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::hes1modelODE(parameters, t(state))))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)

xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- xtrue

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])  
}

xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)

xsim <- insertNaN(xsim.obs,config$filllevel)

tvec.full <- xsim$time
tvec.nobs <- xsim.obs$time

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r.nobs <- abs(foo)
r2.nobs <- r.nobs^2
signr.nobs <- -sign(foo)

cursigma <- rep(NA, ncol(xsim)-1)
curphi <- matrix(NA, 2, ncol(xsim)-1)

for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                    r.nobs, config$kernel)$value
  gr <- function(par) -as.vector(phisigllikC( par, data.matrix(xsim.obs[,1+j]), 
                                              r.nobs, config$kernel)$grad)
  marlikmap <- optim(rep(100, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  cursigma[j] <- marlikmap$par[3]
  curphi[,j] <- marlikmap$par[1:2]
}
cursigma
curphi

j <- 1
# plot(xsim.obs[, "time"], data.matrix(xsim.obs[,1+j]))
phisigllikC(c(curphi[,j], cursigma[j]), data.matrix(xsim.obs[,1+j]), 
            r.nobs, config$kernel)
phisigllikC(c(100, 100, cursigma[j]), data.matrix(xsim.obs[,1+j]), 
            r.nobs, config$kernel)


xtrueAtDiscretization <- sapply(xtrueFunc, function(f) f(xsim[,"time"]))
dotxtrueAtDiscretization = gpds:::hes1modelODE(pram.true$theta, xtrueAtDiscretization)
dotxNumerical <- sapply(1:3, function(j){
  x <- xtrueAtDiscretization[,j]
  dx <- (x[-1] - x[-length(x)])/mean(diff(xsim[,"time"]))
  dx <- apply(cbind(c(NA,dx), c(dx,NA)), 1, mean, na.rm=TRUE)
  dx
})

curphiWithTruth <- matrix(NA, nrow = 2, ncol = 3)
for(j in 1:(ncol(xsim)-1)){
  fn <- function(par) -phillikwithxdotx( par, xtrueAtDiscretization[,j], dotxtrueAtDiscretization[,j],
                                         r, signr, config$kernel)
  marlikmap <- optim(rep(1, 2), fn, method="L-BFGS-B", lower = 0.0001,
                     upper = c(Inf, 60*4*2, Inf))
  
  curphiWithTruth[,j] <- marlikmap$par[1:2]
}
curphiWithTruth  # doesn't work because over-smooth, phi2 is too large

j <- 2
curphi[,j]

fn <- function(par) -phillikwithxdotx( par, xtrueAtDiscretization[,j], dotxNumerical[,j],
                                       r, signr, config$kernel)
marlikmap <- optim(rep(1, 2), fn, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2, Inf))
marlikmap$par  # using dotxNumerical is ok but doesn't use ode information


fn <- function(par) -phillikwithxdotx( par, xtrueAtDiscretization[,j], dotxtrueAtDiscretization[,j],
                                       r, signr, config$kernel)
marlikmap <- optim(rep(100, 2), fn, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2, Inf))
marlikmap$par  # number is way too large, supplying true value makes the system over-smooth


fn <- function(par) -phisigllikC( par, as.matrix(xtrueAtDiscretization[,j]), 
                                  r, config$kernel)$value
gr <- function(par) -as.vector(phisigllikC( par, as.matrix(xtrueAtDiscretization[,j]), 
                                            r, config$kernel)$grad)
marlikmap <- optim(rep(1, 3), fn, gr, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2, Inf))
marlikmap$par

load("hes1image-to-fit-phi.rda")
gpode$xsampled
gpode$fode

j <- 3

fn <- function(par) -phillikwithxdotx( par, colMeans(gpode$xsampled[,,j]), colMeans(gpode$fode[,,j]),
                                       r, signr, config$kernel)
marlikmap <- optim(rep(10, 2), fn, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2, Inf), hessian = TRUE)
marlikmap$par  # doesn't work because over-smooth, same as for curphiWithTruth


fn <- function(par) -phillikwithxdotx( par, gpode$xsampled[2500,,j], gpode$fode[2500,,j],
                                       r, signr, config$kernel)
marlikmap <- optim(rep(10, 2), fn, method="L-BFGS-B", lower = 0.0001,
                   upper = c(Inf, 60*4*2, Inf), hessian = TRUE)
marlikmap$par  # better than curphiWithTruth, but still too large for no reason
fn(c(10, 10))
