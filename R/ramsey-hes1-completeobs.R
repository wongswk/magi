source("R/ramsey-getx0.R")
# First we will define some variable and parameter names

Hes1varnames <- c("P", "M", "H")
Hes1parnames <- c(paste0("theta",1:7))

# and initial conditions and parameters

x0 = log(c(1.438575, 2.037488, 17.90385))
names(x0) <- Hes1varnames

Hes1pars <-  c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
names(Hes1pars) <- Hes1parnames


# The following is a function specifying the model in a form that
# lsoda works with

hes1.ode <- function(times, x, p) {
  dx <- x
  dimnames(dx) <- dimnames(x)
  
  dx["P"] <- -p["theta1"] * exp(x["H"]) + p["theta2"] * exp(x["M"])/exp(x["P"]) - p["theta3"] 
  dx["M"] <- -p["theta4"]  + p["theta5"]/(1 + exp(x["P"])^2)/exp(x["M"])
  dx["H"] <- -p["theta1"] * exp(x["P"]) + p["theta6"]/(1 + exp(x["P"])^2)/exp(x["H"]) - p["theta7"]
  
  return(list(dx))
}

# We need the times at which we will solve the system

times <- seq(0, 60*4, by = 0.01)

# And now we can create solutions to the equations

out <- lsoda(x0, times = times, hes1.ode, Hes1pars)

par(mar = c(5, 5, 1, 1))
matplot(out[, 1], exp(out[, 2:4]), type = "l", xlab = "time", ylab = "(P,M,H)", lwd = 3, 
        cex.lab = 2.5, cex.axis = 2.5)
legend("bottomleft", c("P", "M", "H"), lwd = 3, col = 1:2, lty = 1:2, cex = 1.5)

# We now add some noise to the values of the curves at a reduced set of observation
# times:

hes1times <- seq(0, 240, by=7.5)
nobs <- length(hes1times)
out <- lsoda(x0, times = hes1times, hes1.ode, Hes1pars)
hes1compdata <- out[, 2:4] + 0.15 * matrix(rnorm(3 * nobs), nobs, 3)

# In order to run the profiling proceedures, we need to define some objects.

# The following code will define a basis

hes1range <- c(0, 240)
breaks <- seq(0, 240, 3.75)  ### 3.75 is one break between each time point, for 64 total breaks.
norder <- 4
hes1basis <- create.bspline.basis(range = hes1range, norder = norder, breaks = breaks)

# And this will smooth the data
hes1fdPar <- fdPar(hes1basis, int2Lfd(2), 1)

DEfd0 <- smooth.basis(hes1times, hes1compdata, hes1fdPar)$fd

# and produce a plot of how well this agrees with the smooth
par(mfrow = c(2, 1), mar = c(5, 5, 2, 1))
plotfit.fd(hes1compdata, hes1times, DEfd0, cex.axis = 2.5, cex.lab = 2.5, lwd = 2, cex = 1.5)


# We can now extract the coefficients, which we will also require

coefs0 <- DEfd0$coef
colnames(coefs0) <- Hes1varnames

# CollocInfer requires the right hand side function to be defined somewhat
# differently to lsoda. Here we allow a matrix of values as input

hes1.fun <- function(times, x, p, more) {
  dx <- x
  dx[,"P"] <- -p["theta1"] * exp(x[,"H"]) + p["theta2"] * exp(x[,"M"])/exp(x[,"P"]) - p["theta3"] 
  dx[,"M"] <- -p["theta4"]  + p["theta5"]/(1 + exp(x[,"P"])^2)/exp(x[,"M"])
  dx[,"H"] <- -p["theta1"] * exp(x[,"P"]) + p["theta6"]/(1 + exp(x[,"P"])^2)/exp(x[,"H"]) - p["theta7"]
  
  return(dx)
}


# Now we can choose a trade-off parameter and set up the objects that the
# profiling functions will use.

lambda <- 1000

profile.obj <- LS.setup(pars = Hes1pars, fn = hes1.fun, lambda = lambda, times = hes1times, 
                        coefs = coefs0, basisvals = hes1basis)

proc <- profile.obj$proc
lik <- profile.obj$lik

## Gradient matching can be obtained thr ParsMatchOpt and produces useful initial
## parameter estimates

Pres0 <- ParsMatchOpt(Hes1pars, coefs0, proc)
pars1 <- Pres0$pars  ### rough guess for initial parameters  (based on complete data)

# Profile.LS will do both setup and profiling
# Ores2.2 <- Profile.LS(hes1.fun, hes1compdata, hes1times, Hes1pars, coefs0, hes1basis, lambda)


##### Read incomplete dataset
#data2 <-  read.csv("github/dynamic-systems/R/xsim-obs-Hes1-log.csv")
#data2 <-  as.matrix(data2[,3:5])
#colnames(data2) <- Hes1varnames
data2 <- hes1compdata

### reset smoothing for missing component
#coefs0.2 <- coefs0
#coefs0.2[, 3] <- 0

# The function FitMatchOpt allows us to pull some columns of the coefficient
# matrix into line with the differential equation (keeping the other columns
# fixed):
#Fres3 <- FitMatchOpt(coefs0.2, 3, pars1, proc)  # col 3 is hidden/missing state

### Now run profiling with default lambda = 1000 
#Ores <- Profile.LS(hes1.fun, data2, hes1times, pars1, Fres3$coefs, hes1basis, lambda)
Ores <- Profile.LS(hes1.fun, data2, hes1times, pars1, coefs0, hes1basis, lambda)
def.pars <- Ores$pars  ### parameter estimate.  (with fixed lambda, which might not be ideal)

t.fd <- fd(Ores$coefs, hes1basis)
def.x0 <- eval.fd(0,t.fd)  ## quick hack to get initial B-spline-fitted x0

# Covariance of parameters can be obtained from
covar <- Profile.covariance(Ores$pars, times = hes1times, data = data2, coefs = Ores$coefs, 
                            lik = lik, proc = proc)
covar

# and we can look at confidence intervals
CIs <- cbind(Ores$pars - 2 * sqrt(diag(covar)), Ores$pars + 2 * sqrt(diag(covar)))
rownames(CIs) <- Hes1parnames
CIs


#### We can run with lambda ladder to find optimal tuning parameter

whichtimes <- cbind(1:25, 9:33)  ### break up data for tuning lambda based on forward prediction error

FPE <- forward.prediction.error(hes1times, data2, Ores$coefs, lik, 
                                proc, Ores$pars, whichtimes)  # based on default lambda = 1000


par(mar = c(5, 5, 1, 1))
matplot(hes1times, data2, pch = Hes1varnames, cex = 1.5, cex.lab = 2.5, cex.axis = 2.5)

lambdas <- 10^(3:7)
FPEs <- 0 * lambdas
all.pars <- matrix(NA, nrow = length(lambdas), ncol = length(pars1))
all.x0 <- matrix(NA, nrow = length(lambdas), ncol = length(x0))
temp.pars <- pars1
temp.coefs <- Fres3$coefs
for (ilam in 1:length(lambdas)) {
  print(paste("lambda = ", lambdas[ilam]))
  t.Ores <- Profile.LS(hes1.fun, data2, hes1times, temp.pars, temp.coefs, hes1basis, lambdas[ilam])
  print(t.Ores$pars)
  temp.pars <- t.Ores$pars
  all.pars[ilam,] <- t.Ores$pars
  temp.coefs <- t.Ores$coefs
  
  t.fd <- fd(t.Ores$coefs, hes1basis)
  lines(t.fd, lwd = 2, col = ilam, lty = 1)
  
  all.x0[ilam,] <- eval.fd(0, t.fd)  ## quick hack to get initial x0
  
  FPEs[ilam] <- forward.prediction.error(hes1times, data2, t.Ores$coefs, lik, 
                                         proc, t.Ores$pars, whichtimes)
}

#legend("bottomleft", c("lambda = 1", "lambda = 1e3", "lambda = 1e7"), col = c(1,  2, 4), lwd = 2, cex = 1.2)

FPEs

par(mar = c(5, 5, 1, 1))
plot(log10(lambdas), FPEs, type = "l", lwd = 2, col = 4, cex.lab = 2.5, cex.axis = 2.5)

best.pars <- all.pars[which.min(FPEs),]   ## can pick lowest FPE as the parameter estimate and associated x0
best.pars.x0 <- all.x0[which.min(FPEs),]

