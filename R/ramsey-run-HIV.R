# First we will define some variable and parameter names
library(deSolve)
library(CollocInfer)


#### all possible data files

# rdaDir <- "../results/cpp/fullobs/"
 
subdirs <- c("../comparison/results/HIV-9obs-fill3")
#              "../results/cpp/fullobs//variablephi-temper-coldstart", "../results/cpp/fullobs//variablephi-temper-warmstart")
all_files <- lapply(subdirs, list.files)
all_files <- lapply(all_files, function(x) x[grep(".*HIV-([0-9]+)-noise.*rda", x)])

all_seeds <- lapply(all_files, function(x) gsub(".*HIV-([0-9]+)-noise.*rda", "\\1", x))

all_seeds <- all_seeds[[1]]
# common_seeds <- all_seeds[[1]]
# for (i in 2:length(all_seeds)){
#   common_seeds <- intersect(common_seeds, all_seeds[[i]])  
# }

dir.create(paste0(subdirs[1], "/ramsay"), showWarnings = FALSE)



# for (f in rda_files){
#   if (!file.exists(paste0(rdaDir, "ramsay-optimation-", f))){
#     print(f)
#   }
# }

# args <- commandArgs(trailingOnly = TRUE)
# args <- as.numeric(args)
# if(length(args) == 0){
#   rda_it = 1
# }else{
#   rda_it = args
# }

for (rda_it in 1:length(all_seeds)) {
#rda_it <- 1

envHIV <- new.env()
##### Read dataset
each_seed <- common_seeds[rda_it]

# if (file.exists(paste0("../results/cpp/fullobs/ramsay/Hes1-log-",each_seed,"-fullobs-ramsay.rda"))){
#   quit(save = "no")
# }
# 


# system(paste0(" rm -r ../fullobs", each_seed))
# dir.create(paste0("../fullobs", each_seed), showWarnings = FALSE)
# setwd(paste0("../fullobs", each_seed))
# 
# pdf(paste0("../results/cpp/fullobs/ramsay/Hes1-log-",each_seed,"-fullobs-ramsay.pdf"))
# sink(paste0("../results/cpp/fullobs/ramsay/Hes1-log-",each_seed,"-fullobs-ramsay.txt"))

print(each_seed)
load(paste0(subdirs[1], "/HIV-",each_seed,"-noise0.3.rda"), envir = envHIV)

data2 <-  envHIV$xsim[complete.cases(envHIV$xsim),]
data2 <-  as.matrix(data2[,-1])


### smoothing based on the observed data -----
HIVvarnames <- c("T", "Tm", "Tw", "Tmw")
HIVparnames <- c(paste0("theta",1:9))
colnames(data2) <- HIVvarnames

# and initial conditions and parameters

x0 = log(c(3.35e7, 134000, 25000, 9000))
names(x0) <- HIVvarnames

HIVpars <-  c(0.015, 1.51e-3, 1.11e-3, 4.4e-4, 4.15e-3, 1.1e-3, -2.29e-2, 7.13e-03, 5.68e-04)
names(HIVpars) <- HIVparnames


# The following is a function specifying the model in a form that
# lsoda works with

HIV.ode <- function(times, x, p) {
  dx <- x
  dimnames(dx) <- dimnames(x)

  T = exp(x["T"])
  Tm = exp(x["Tm"])
  Tw = exp(x["Tw"])
  Tmw = exp(x["Tmw"])
  
  dx["T"] <- (p["theta1"] - 1e-6*p["theta2"]*Tm - 1e-6*p["theta3"]*Tw - 1e-6*p["theta4"]*Tmw)
  dx["Tm"] <- (p["theta7"] + 1e-6*p["theta2"]*T - 1e-6*p["theta5"]*Tw) + 1e-6*0.25*p["theta4"]*Tmw*T / Tm
  dx["Tw"] <- (p["theta8"] + 1e-6*p["theta3"]*T - 1e-6*p["theta6"]*Tm) + 1e-6*0.25*p["theta4"]*Tmw*T / Tw
  dx["Tmw"] <- p["theta9"] + 0.5*1e-6*p["theta4"]*T + (1e-6*p["theta5"]+1e-6*p["theta6"])*Tw*Tm / Tmw

  return(list(dx))
}

# We need the times at which we will solve the system

times <- seq(0, 93, by = 0.25)

# And now we can create solutions to the equations

out <- lsoda(x0, times = times, HIV.ode, HIVpars)

par(mar = c(5, 5, 1, 1))
matplot(out[, 1], exp(out[, 2:5]), type = "l", xlab = "time", ylab = "(T,Tm,Tw,Tmw)", lwd = 3, 
        cex.lab = 2.5, cex.axis = 2.5)
legend("bottomleft", HIVvarnames, lwd = 3, col = 1:2, lty = 1:2, cex = 1.5)

# We now add some noise to the values of the curves at a reduced set of observation
# times:

HIVtimes <- c(70,82, 94,106,115,127,139,151,163)-70
nobs <- length(HIVtimes)
out <- lsoda(x0, times = HIVtimes, HIV.ode, HIVpars)
HIVdata <- data2  ## read from rda so noise already added for that seed

# In order to run the profiling proceedures, we need to define some objects.

# The following code will define a basis

HIVrange <- c(0, 93)
breaks <- seq(0, 93, 3)  ### 3.75 is one break between each time point, for 64 total breaks.
norder <- 4
HIVbasis <- create.bspline.basis(range = HIVrange, norder = norder, breaks = breaks)

# And this will smooth the data
HIVfdPar <- fdPar(HIVbasis, int2Lfd(2), 1)

DEfd0 <- smooth.basis(HIVtimes, HIVdata, HIVfdPar)$fd

# and produce a plot of how well this agrees with the smooth
# par(mfrow = c(2, 1), mar = c(5, 5, 2, 1))
# plotfit.fd(HIVdata, HIVtimes, DEfd0, cex.axis = 2.5, cex.lab = 2.5, lwd = 2, cex = 1.5)


# We can now extract the coefficients, which we will also require

coefs0 <- DEfd0$coef
colnames(coefs0) <- HIVvarnames

# CollocInfer requires the right hand side function to be defined somewhat
# differently to lsoda. Here we allow a matrix of values as input

HIV.fun <- function(times, x, p, more) {
  dx <- x
  
  T = exp(x[,"T"])
  Tm = exp(x[,"Tm"])
  Tw = exp(x[,"Tw"])
  Tmw = exp(x[,"Tmw"])
  
  dx[,"T"] <- (p["theta1"] - 1e-6*p["theta2"]*Tm - 1e-6*p["theta3"]*Tw - 1e-6*p["theta4"]*Tmw)
  dx[,"Tm"] <- (p["theta7"] + 1e-6*p["theta2"]*T - 1e-6*p["theta5"]*Tw) + 1e-6*0.25*p["theta4"]*Tmw*T / Tm
  dx[,"Tw"] <- (p["theta8"] + 1e-6*p["theta3"]*T - 1e-6*p["theta6"]*Tm) + 1e-6*0.25*p["theta4"]*Tmw*T / Tw
  dx[,"Tmw"] <- p["theta9"] + 0.5*1e-6*p["theta4"]*T + (1e-6*p["theta5"]+1e-6*p["theta6"])*Tw*Tm / Tmw
  
  return(dx)
}

# Now we can choose a trade-off parameter and set up the objects that the
# profiling functions will use.

lambda <- 1000

profile.obj <- LS.setup(pars = HIVpars, fn = HIV.fun, lambda = lambda, times = HIVtimes, 
                        coefs = coefs0, basisvals = HIVbasis)

proc <- profile.obj$proc
lik <- profile.obj$lik

## Gradient matching can be obtained thr ParsMatchOpt and produces useful initial
## parameter estimates

Pres0 <- ParsMatchOpt(HIVpars, coefs0, proc)
pars1 <- Pres0$pars  ### rough guess for initial parameters  (based on complete data)

# Profile.LS will do both setup and profiling
# Ores2.2 <- Profile.LS(hes1.fun, hes1compdata, hes1times, Hes1pars, coefs0, hes1basis, lambda)

### Now run profiling with default lambda = 1000 
Ores <- Profile.LS(HIV.fun, data2, HIVtimes, pars1, coefs0, HIVbasis, lambda)
#Ores$pars <- pmax(Ores$pars, 1e-3)  # need positive parameters to run numerical solver
def.pars <- Ores$pars  ### parameter estimate.  (with fixed lambda, which might not be ideal)

t.fd <- fd(Ores$coefs, HIVbasis)
#def.x0 <- getfdx0(t.fd)  ## quick hack to get initial B-spline-fitted x0
def.x0 <- eval.fd(0, t.fd)

# Covariance of parameters can be obtained from
covar <- Profile.covariance(Ores$pars, times = HIVtimes, data = data2, coefs = Ores$coefs, 
                            lik = lik, proc = proc)
covar

# and we can look at confidence intervals
CIs <- cbind(Ores$pars - 2 * sqrt(diag(covar)), Ores$pars + 2 * sqrt(diag(covar)))
rownames(CIs) <- HIVparnames
CIs


#### We can run with lambda ladder to find optimal tuning parameter

whichtimes <- cbind(1:6, 4:9)  ### break up data for tuning lambda based on forward prediction error

FPE <- forward.prediction.error(HIVtimes, data2, Ores$coefs, lik, 
                                proc, Ores$pars, whichtimes)  # based on default lambda = 1000


par(mar = c(5, 5, 1, 1))
matplot(HIVtimes, data2, pch = HIVvarnames, cex = 1.5, cex.lab = 2.5, cex.axis = 2.5)

lambdas <- 10^(3:7)
FPEs <- 0 * lambdas
all.pars <- matrix(NA, nrow = length(lambdas), ncol = length(pars1))
all.x0 <- matrix(NA, nrow = length(lambdas), ncol = length(x0))
all.Ores <- list()
all.ci <- list()
temp.pars <- pars1
temp.coefs <- coefs0
for (ilam in 1:length(lambdas)) {
  print(paste("lambda = ", lambdas[ilam]))
  t.Ores <- Profile.LS(HIV.fun, data2, HIVtimes, temp.pars, temp.coefs, HIVbasis, lambdas[ilam])
  #t.Ores$pars <- pmax(t.Ores$pars, 1e-3)  # need positive parameters to run numerical solver
  all.Ores[[ilam]] <- t.Ores
  print(t.Ores$pars)
  temp.pars <- t.Ores$pars
  all.pars[ilam,] <- t.Ores$pars
  temp.coefs <- t.Ores$coefs
  
  t.fd <- fd(t.Ores$coefs, HIVbasis)
  lines(t.fd, lwd = 2, col = ilam, lty = 1)
  
  #all.x0[ilam,] <- getfdx0(t.fd)  ## quick hack to get initial x0
  all.x0[ilam,] <- eval.fd(0, t.fd)
  
  FPEs[ilam] <- forward.prediction.error(HIVtimes, data2, t.Ores$coefs, lik, 
                                         proc, t.Ores$pars, whichtimes)
  
  # Covariance of parameters can be obtained from
  covar <- Profile.covariance(t.Ores$pars, times = HIVtimes, data = data2, coefs = t.Ores$coefs, 
                              lik = lik, proc = proc)
  covar
  
  # and we can look at confidence intervals
  CIs <- cbind(t.Ores$pars - 2 * sqrt(diag(covar)), t.Ores$pars + 2 * sqrt(diag(covar)))
  rownames(CIs) <- HIVparnames
  all.ci[[ilam]] <- CIs
}

#legend("bottomleft", c("lambda = 1", "lambda = 1e3", "lambda = 1e7"), col = c(1,  2, 4), lwd = 2, cex = 1.2)

FPEs

par(mar = c(5, 5, 1, 1))
plot(log10(lambdas), FPEs, type = "l", lwd = 2, col = 4, cex.lab = 2.5, cex.axis = 2.5)

best.pars <- all.pars[which.min(FPEs),]   ## can pick lowest FPE as the parameter estimate and associated x0
best.pars.x0 <- all.x0[which.min(FPEs),]
best.ci <- all.ci[[which.min(FPEs)]]

modelODE <- function(t, state, parameters) {
  list(as.vector(gpds:::HIVmodelODE(parameters, t(state))))
}

xdesolveRamsay <- deSolve::ode(y = best.pars.x0, times = times, func = modelODE, parms = best.pars)
# matplot(xdesolveRamsay[,1], xdesolveRamsay[,-1], col=1:3, type="l")

rmse_log <- sqrt(colMeans((out[,-1] - xdesolveRamsay[match(out[,1], xdesolveRamsay[,1]),-1])^2))
#rmse_orig <- sqrt(colMeans((exp(out[,-1]) - exp(xdesolveRamsay[match(out[,1], xdesolveRamsay[,1]),-1]))^2))

save(list=setdiff(ls(), "envHIV"), file = paste0(subdirs[1], "/ramsay/HIV-",each_seed,"-ramsay.rda"))
#dev.off()
#setwd("../dynamic-systems/")

}

