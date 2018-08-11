library(gpds)
library(rstan)
options(mc.cores = parallel::detectCores())

# simulation set up ------------------------------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 101,
    noise = c(0.001, 0.001),
    seed = 125455454, #(as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,
    npostplot = 50,
    filllevel = 1,
    modelName = "FN",
    kernel = "generalMatern-4.5"
  )
}

if(!exists("stanConfig")){
  stanConfig <- list(
    sigma_obs=0.5,
    sigma_xdot=0.1
  )
}

config$ndis <- (config$nobs-1)*2^config$filllevel+1

if(grepl("/n/",getwd())){
  baseDir <- "/n/regal/kou_lab/shihaoyang/DynamicSys/results/" # tmp folder on cluster 
  # config$seed <- (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9 # random seed on cluster
}else{
  baseDir <- "~/Workspace/DynamicSys/results/batch-output/"  
}

outDir <- with(config, paste0(baseDir,"stan/"))
system(paste("mkdir -p", outDir))
filename <- paste0(outDir, 
                   "nobs-",config$nobs,"_",
                   "filllevel-",config$filllevel,"_",
                   config$kernel,"_",
                   "noise-",paste(config$noise, collapse = "-"),"_",
                   "seed",config$seed)

pram.true <- list(
  theta=c(0.2,0.2,3),
  x0 = c(-1, 1),
  phi=c(0.9486433, 3.2682434,
        1.9840824, 1.1185157),
  sigma=config$noise
)

times <- seq(0,20,0.05)

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::fnmodelODE(parameters, t(state))))
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

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

xsim <- insertNaN(xsim.obs,config$filllevel)
obs_index <- match(rownames(xsim.obs), rownames(xsim))

xtrue.atsim <- sapply(xtrueFunc, function(f) f(xsim$time))
dotxtrue.atsim <- gpds:::fnmodelODE(pram.true$theta, xtrue.atsim)

dotxtrue = gpds:::fnmodelODE(pram.true$theta, data.matrix(xtrue[,-1]))

init <- list(
  abc = pram.true$theta,
  vtrue = xtrue.atsim[,1],
  rtrue = xtrue.atsim[,2]
)

config$stanConfig <- stanConfig

curCov <- getCovMphi(config$kernel, xsim, xsim.obs)

pdf(paste0(filename, "_m-visualize.pdf"))
image(curCov[[1]]$mphi, col = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(n = 101),
      main="m-matrix for component V")
image(curCov[[1]]$mphi, col = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(n = 101),
      main="m-matrix for component R")
dev.off()

# STAN sampling using finite difference ----------------------------------------

gpsmooth <- stan(file="stan/m-only.stan",
                 data=list(
                   n_discret=nrow(xsim),
                   n_obs=nrow(xsim.obs),
                   vobs=xsim.obs$X1,
                   robs=xsim.obs$X2,
                   time=xsim$time,
                   obs_index=obs_index,
                   sigma_obs=stanConfig$sigma_obs,
                   sigma_xdot=stanConfig$sigma_xdot,
                   mphi_v=curCov[[1]]$mphi,
                   mphi_r=curCov[[2]]$mphi
                 ),
                 iter=1000, chains=7, init=rep(list(init), 7), warmup = 300, 
                 control = list(max_treedepth = 15))

gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

stanode <- list(theta=gpsmooth_ss$abc,
                xsampled=abind::abind(gpsmooth_ss$vtrue, gpsmooth_ss$rtrue, along = 3),
                lglik=gpsmooth_ss$lp__,
                sigma = pram.true$sigma)
stanode$fode <- sapply(1:length(stanode$lglik), function(t) 
  with(stanode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
stanode$fode <- aperm(stanode$fode, c(3,1,2))


gpds:::plotPostSamplesFlex(
  paste0(filename, "_m-stan.pdf"), 
  xtrue, dotxtrue, xsim, stanode, pram.true, config)

save.image(paste0(filename, "_m-stan.rda"))

# likelihood move away from mode ------------------------------------------------
gpsmooth <- stan(file="stan/m-only.stan",
                 data=list(
                   n_discret=nrow(xsim),
                   n_obs=nrow(xsim.obs),
                   vobs=xsim.obs$X1,
                   robs=xsim.obs$X2,
                   time=xsim$time,
                   obs_index=obs_index,
                   sigma_obs=stanConfig$sigma_obs,
                   sigma_xdot=stanConfig$sigma_xdot,
                   mphi_v=curCov[[1]]$mphi,
                   mphi_r=curCov[[2]]$mphi
                 ),
                 iter=100, chains=1, init=rep(list(init), 1), warmup = 0)

gpsmooth_ss <- extract(gpsmooth, permuted=TRUE)

stanode <- list(theta=gpsmooth_ss$abc,
                xsampled=abind::abind(gpsmooth_ss$vtrue, gpsmooth_ss$rtrue, along = 3),
                lglik=gpsmooth_ss$lp__,
                sigma = pram.true$sigma)
stanode$fode <- sapply(1:length(stanode$lglik), function(t) 
  with(stanode, gpds:::fnmodelODE(theta[t,], xsampled[t,,])), simplify = "array")
stanode$fode <- aperm(stanode$fode, c(3,1,2))

gpds:::plotPostSamplesFlex(
  paste0(filename, "_m-likelihoodmove-stan.pdf"), 
  xtrue, dotxtrue, xsim, stanode, pram.true, config)

system(paste0("pdftk ",
              filename, "_m-likelihoodmove-stan.pdf ",
              " cat 8 output ",
              filename, "_m-likelihoodtrue-stan.pdf "
))

system(paste0("pdftk ", filename, "_m-stan.pdf ", 
              filename, "_m-likelihoodtrue-stan.pdf ",
              filename, "_m-visualize.pdf ",
              " cat output ",
              filename, "_m-comb.pdf"
              ))
