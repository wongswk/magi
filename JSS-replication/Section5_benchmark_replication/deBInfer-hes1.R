# Run deBInfer method on Hes1 model

args <- commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
if(length(args) == 0){
  seed = 1
}else{
  seed = args
}

dir.create("results/deBInfer", showWarnings = FALSE)
pdf(paste0("results/deBInfer/Hes1-log-",seed,".pdf"))
sink(paste0("results/deBInfer/Hes1-log-",seed,".txt"))

print(seed)
set.seed(seed)

source("setup-sim.R")

library(deBInfer)

StartTime <- proc.time()[3]

modelR <- function(t, Y, parameters) {
  with (as.list(parameters),{
    dP <- -par_a * exp(Y[3]) + par_b * exp(Y[2])/ exp(Y[1]) - par_c
    dM <- -par_d + par_e / (1 + exp(Y[1])^2) / exp(Y[2])
    dH <- -par_a * exp(Y[1]) + par_f/(1 + exp(Y[1])^2)/exp(Y[3]) - par_g
    list(c(dP, dM, dH))
  })
}

jacR <- function (t, Y, parameters) {
  with (as.list(parameters),{
    dP = -(1 + exp(2 * Y[1]))^(-2) * exp(2 * Y[1]) * 2
    
    PD[1,1] <- -par_b * exp(Y[2]-Y[1])
    PD[1,2] <- par_b * exp(Y[2]-Y[1])
    PD[1,3] <- -par_a * exp(Y[3])
    PD[2,1] <- par_e * exp(-Y[2]) * dP
    PD[2,2] <- -par_e * exp(-Y[2]) / (1 + exp(2 * Y[1]))
    PD[2,3] <- 0
    PD[3,1] <- -par_a * exp(Y[1]) + par_f * exp(-Y[3]) * dP
    PD[3,2] <- 0
    PD[3,3] <- -par_f * exp(-Y[3]) / (1 + exp(2 * Y[1]))
    return(PD)
  })
}

# the observation model
obs_model <- function(data, sim.data, samp){
  P_ind <- y$time %in% seq(0, 240, by = 15)
  M_ind <- y$time %in% seq(7.5, 240, by = 15)
  
  llik.P <- sum(dnorm(y[P_ind,2], mean = sim.data[P_ind,2], sd = samp[['sd.P']], log = TRUE))
  llik.M <- sum(dnorm(y[M_ind,3], mean = sim.data[M_ind,3], sd = samp[['sd.M']], log = TRUE))
  return(llik.P + llik.M)
}


par_a <- debinfer_par(name = "par_a", var.type = "de", fixed = FALSE,
                   value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 2),
                   prop.var = 0.01, samp.type = "rw-ref")
par_b <- debinfer_par(name = "par_b", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 2),
                      prop.var = 0.01, samp.type = "rw-ref")
par_c <- debinfer_par(name = "par_c", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 2),
                      prop.var = 0.01, samp.type = "rw-ref")
par_d <- debinfer_par(name = "par_d", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 2),
                      prop.var = 0.01, samp.type = "rw-ref")
par_e <- debinfer_par(name = "par_e", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 2),
                      prop.var = 0.05, samp.type = "rw-ref")
par_f <- debinfer_par(name = "par_f", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 50), prior = "unif", hypers = list(min = 0, max = 100),
                      prop.var = 0.5, samp.type = "rw-ref")
par_g <- debinfer_par(name = "par_g", var.type = "de", fixed = FALSE,
                      value = runif(1, 0, 1), prior = "unif", hypers = list(min = 0, max = 10),
                      prop.var = 0.01, samp.type = "rw-ref")

sd.P <- debinfer_par(name = "sd.P", var.type = "obs", fixed = TRUE, value = 0.15)
sd.M <- debinfer_par(name = "sd.M", var.type = "obs", fixed = TRUE, value = 0.15)


P <- debinfer_par(name = "P", var.type = "init", fixed = FALSE,
                  value = log(param.true$x0[1]), prior = "unif", hypers = list(min = -10, max = 10),
                  prop.var = 0.05, samp.type = "rw-ref")
M <- debinfer_par(name = "M", var.type = "init", fixed = FALSE,
                  value = log(param.true$x0[2]), prior = "unif", hypers = list(min = -10, max = 10),
                  prop.var = 0.05, samp.type = "rw-ref")
H <- debinfer_par(name = "H", var.type = "init", fixed = FALSE,
                  value = log(param.true$x0[3]), prior = "unif", hypers = list(min = -10, max = 10),
                  prop.var = 0.05, samp.type = "rw-ref")

mcmc.pars <- setup_debinfer(par_a, par_b, par_c, par_d, par_e, par_f, par_g, sd.P, sd.M, P, M, H)

iter <- 20000
mcmc_samples <- de_mcmc(N = iter, data = y, de.model = modelR,
                        obs.model = obs_model, all.params = mcmc.pars,
                        Tmax = max(y[,"time"]), data.times = y[,"time"], cnt = iter*0.05,
                        plot = TRUE, solver = "ode", verbose.mcmc = TRUE)                        

TimeUsed <- proc.time()[3] - StartTime


burnin <- 0.5  # burn-in proportion

tmp_samples <- mcmc_samples$samples[ (burnin*iter+1):iter,]
tmp_lp <- mcmc_samples$lpost[ (burnin*iter+1):iter]

theta.est <- apply(tmp_samples[,1:7], 2, mean)
x0.est <- apply(tmp_samples[,8:10], 2, mean)

#### Calculate trajectory RMSE at the 33 time points
xdesolve_recon <- deSolve::ode(y = x0.est, times = seq(0, 240, by = 7.5), func = logmodelODE, parms = theta.est)

rmse_log <- sqrt(colMeans((x[match(xdesolve_recon[,1], x[,1]),-1] - xdesolve_recon[,-1])^2))
rmse_orig <- sqrt(colMeans((exp(x[match(xdesolve_recon[,1], x[,1]),-1]) - exp(xdesolve_recon[,-1]))^2))

print(rmse_log)
print(rmse_orig)

save.image(file = paste0("results/deBInfer/Hes1-log-",seed,".rda"))

