library(magi)

outDir <- "../results/Michaelis-Menten-Va/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
realdata <- read.csv(paste0("../results/Michaelis-Menten/", "hydrolysis.csv"))
realdata$P = realdata$P / 2

noise_scalar = 0.02
seed = 1

noise = c(NA, noise_scalar, noise_scalar)
linfillspace = c(0.5)
linfillcut = NULL
phi = cbind(c(0.1, 70), c(1, 30), c(1, 30))
phi_change_time = 0
time_acce_factor = 1
obs_keep = 1:26
obs_source = "va-csv"
t.truncate = 70



# set up configuration if not already exist ------------------------------------

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
  obs_source = obs_source,
  modelName = "Michaelis-Menten-Va"
)

if(is.null(config$skip_visualization)){
  matplot(realdata$t, realdata[,-1], type="b")
}

# initialize global parameters, true x, simulated x ----------------------------
# parameters and initial conditions that seem to mimic the real data well
pram.true <- list(
  theta=c(0.35, 0.2, 2.54),
  # x0 = c(0.08277011, 0.7500533, 0.2327168),
  x0 = c(0.1, 1, 0),
  phi = config$phi,
  sigma=config$noise
)

times <- seq(0,config$t.end,0.01)

modelODE <- function(t, state, parameters) {
  list(as.vector(magi:::MichaelisMentenModelVaODE(parameters, t(state), t)))
}

get_sse <- function(theta){
  xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = theta)
  xtrue <- data.frame(xtrue)
  idx = match(round(realdata$t, 2), round(xtrue$time, 2))
  sum((xtrue[idx,c(4,3)] - realdata[,-1])^2)
}
get_sse(pram.true$theta)
get_sse(c(1.170971, 0, 0.324367))
get_sse(c(0.2752736, -0.2956150,  0.4217417))
get_sse(c(1.3855510, 0.1000000, 0.3578312))

optim_result <- optim(pram.true$theta, get_sse, lower=0.1)


theta = c(1.17, 0, 0.32)

xtrue <- deSolve::ode(y = pram.true$x0, times = times, func = modelODE, parms = theta)
xtrue <- data.frame(xtrue)
idx = match(round(realdata$t, 2), round(xtrue$time, 2))
sum((xtrue[idx,c(4,3)] - realdata[,-1])^2)

if(is.null(config$skip_visualization)){
  matplot(xtrue[, "time"], xtrue[, c(3,4)], type="l", lty=1)
  matplot(realdata$t, realdata[,-1], type="p", add=TRUE)
}


