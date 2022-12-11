library(magi)


outDir <- "../results/Michaelis-Menten-Vb4p/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0("../results/Michaelis-Menten/", "hydrolysis.csv"))
realdata$P = realdata$P / 2

noise_scalar = 0.02
outDir <- paste0("../results/Michaelis-Menten-Vb4p/", "noise", noise_scalar, "/")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

seed = 1

noise = c(NaN, noise_scalar, NaN, noise_scalar)
linfillspace = c(0.5)
linfillcut = NULL
phi = cbind(c(0.1, 70), c(1, 30), c(0.1, 70), c(0.5, 30))
phi_change_time = 0
time_acce_factor = 1
obs_keep = setdiff(1:26, c(1,2,4,6,8,11))
obs_source = "vb-sim"
t.truncate = 70

config <- list(
  nobs = nrow(realdata),
  noise = noise,
  kernel = "generalMatern",
  seed = seed,
  bandsize = 40,
  hmcSteps = 100,
  n.iter = 8001,
  linfillspace = linfillspace, 
  linfillcut = linfillcut,
  t.end = 70,
  t.start = 0,
  obs_start_time = 0,
  phi_change_time = phi_change_time,
  time_acce_factor = time_acce_factor,
  t.truncate = t.truncate,
  obs_keep = obs_keep,
  useMean = TRUE,
  phi = phi,
  skip_visualization = TRUE,
  obs_source = obs_source,
  modelName = "Michaelis-Menten-Vb4p"
)


# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list( 
  theta=c(0.636, 0.0, 10.8, 0.0),
  x0 = c(0.1, 1, 0, 0),
  phi = config$phi
)

times <- seq(0,config$t.end,0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelVb4pODE(parameters, t(state), t)))
}

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = pram.true$theta)
xtrue <- data.frame(xtrue)
sum(log(xtrue[-1,]))
matplot(xtrue[, "time"], xtrue[, -1], type="l", lty=1)
matplot(xtrue[, "time"], xtrue[, c(3,5)], type="l", lty=1)
matplot(realdata$t, realdata[,-1], type="p", add=TRUE)

get_sse <- function(theta){
  xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = theta)
  xtrue <- data.frame(xtrue)
  idx = match(round(realdata$t, 2), round(xtrue$time, 2))
  sum((xtrue[idx,c(5,3)] - realdata[,-1])^2)
}

get_sse(pram.true$theta)

optim_result <- optim(pram.true$theta, get_sse, lower=0)
optim_result$par
optim_result$value
