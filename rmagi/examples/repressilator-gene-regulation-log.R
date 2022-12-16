library(magi)

outDir <- "../results/repressilator-gene-regulation-log/"
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# set up configuration if not already exist ------------------------------------
if(!exists("config")){
  config <- list(
    nobs = 51,
    noise = rep(0.3, 6),
    kernel = "generalMatern",
    seed = 142249801,  # example seed
    #seed = (as.integer(Sys.time())*104729+sample(1e9,1))%%1e9,  # or choose a random seed
    niterHmc = 10001,
    filllevel = 1,
    t.end = 300,
    modelName = "repressilator-gene-regulation-log"
  )
}


# initialize global parameters, true x, simulated x ----------------------------
alpha <- 240 # obtain from Fig 1b in Elowitz and Leibler (2000)
KM <- 40     # scale factor only, to convert protein number to match Fig 1c in paper
pram.true <- list(
  theta=c(0.001*alpha, alpha, 2, 1/5),  # alpha0/alpha = 0.001
  # initial condition cannot be the same, otherwise the system degenerates to two-components 
  # i.e., all the m and all the p will be the same
  x0 = log(c(0.4, 20, 40, 0.01, 0.01, 0.01)),
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
matplot(xtrue[, "time"], exp(xtrue[, -c(1,5,6,7)]), type="l", lty=1)


xtrueFunc <- lapply(2:ncol(xtrue), function(j)
  approxfun(xtrue[, "time"], xtrue[, j]))

xsim <- data.frame(time = seq(0,config$t.end,length=config$nobs))
xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$time)))

set.seed(config$seed)
for(j in 1:(ncol(xsim)-1)){
  xsim[,1+j] <- xsim[,1+j]+rnorm(nrow(xsim), sd=config$noise[j])
}


xsim.obs <- xsim[seq(1,nrow(xsim), length=config$nobs),]
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

# pdf(width = 20, height = 5, file=paste0(outDir, "sample-data.pdf"))
# compnames <- c("m_laci", "m_tetr", "m_ci", "p_laci (unobserved)", "p_tetr (unobserved)", "p_ci (unobserved)")
# par(mfrow=c(1,3))
# for (ii in 1:3) {
#   plot( c(min(xsim.obs$time),max(xsim.obs$time)), c(min(exp(xsim.obs[,ii+1])), max(exp(xsim.obs[,ii+1]))), type='n',xlab="time", ylab='')
#   lines(xtrue[ xtrue$time >=1, "time"], exp(xtrue[xtrue$time >=1,ii+1]),col=ii)
#   mtext(compnames[ii], cex=1.25)
#   points(xsim.obs$time[-1], exp(xsim.obs[-1,ii+1]), col=ii)
# }
# dev.off()

xsim <- setDiscretization(xsim.obs,config$filllevel)

dynamicalModelList <- list(
  fOde=magi:::repressilatorGeneRegulationLogODE,
  fOdeDx=magi:::repressilatorGeneRegulationLogDx,
  fOdeDtheta=magi:::repressilatorGeneRegulationLogDtheta,
  thetaLowerBound=rep(0, 4),
  thetaUpperBound=rep(Inf, 4),
  name="repressilator-gene-regulation-log"
)

testDynamicalModel(dynamicalModelList$fOde, dynamicalModelList$fOdeDx, dynamicalModelList$fOdeDtheta, "dynamicalModelList", data.matrix(xsim.obs[-1,-1]), pram.true$theta, xsim.obs$time[-1])

matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim)-1), pch=20)

# set some reasonable hyperparameters
phiExogenous <- rbind(rep(6, 6), rep(10, 6))
# known noise for mRNA
sigmaInit <- config$noise

# remove initial conditions
xsim <- xsim[-1,]

# set protein levels missing
xsim[,5:7] <- NA
xsim.obs[,5:7] <- NA


OursStartTime <- proc.time()[3] 

result <- magi::MagiSolver(xsim[,-1], dynamicalModelList, xsim$time, 
                           control = list(niterHmc=config$niterHmc, phi=phiExogenous, sigma=sigmaInit, useFixedSigma=TRUE, verbose=TRUE))
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
magi:::plotPostSamplesFlex(
  paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], "-partobs-fill",config$filllevel,".pdf"),
  xtrue, dotxtrue, xsim, gpode, pram.true, config, odemodel)

save.image(paste0(outDir, config$modelName,"-",config$seed,"-noise", config$noise[1], ".rda"))

#### Plot noisy observations and inferred trajectories

## remove initial conditions
xtrue <- xtrue[xtrue$time >= 1,]
xsampledexp <- exp(gpode$xsampled)
oursPostExpX <- cbind(
  apply(xsampledexp, 2:3, mean),
  apply(xsampledexp, 2:3, function(x) quantile(x, 0.025)),
  apply(xsampledexp, 2:3, function(x) quantile(x, 0.975)))

compnames <- c("m_lacI", "m_tetR", "m_cI", expression(paste("p_lacI (", bold("unobserved"), ")")), expression(paste("p_tetR (", bold("unobserved"), ")")), expression(paste("p_cI (", bold("unobserved"), ")")))
layout(rbind(c(1,2,3), c(4,5,6), c(7,7,7)), heights = c(8,8,1))
for (ii in 1:6) {
  
  par(mar = c(4, 4.5, 1.75, 0.1))
  ourEst <- oursPostExpX[,ii]
  ourEst <- exp(magi:::getMeanCurve(xsim$time, log(ourEst), xtrue[,1],
                                    t(phiExogenous[,ii]), 0, 
                                    kerneltype=config$kernel, deriv = FALSE))
  
  ourUB <- oursPostExpX[,12+ii]
  ourUB <- exp(magi:::getMeanCurve(xsim$time, log(ourUB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))
  
  
  ourLB <- oursPostExpX[,6+ii]
  ourLB <- exp(magi:::getMeanCurve(xsim$time, log(ourLB), xtrue[,1],
                                   t(phiExogenous[,ii]), 0,
                                   kerneltype=config$kernel, deriv = FALSE))
  plot( c(min(xtrue$time),max(xtrue$time)), c(min(ourLB), min(max(ourUB),175)), type='n', xlab='', ylab='')
  
  polygon(c(xtrue[,1], rev(xtrue[,1])), c(ourUB, rev(ourLB)),
          col = "skyblue", border = NA)    
  # polygon(c(xsim[,1], rev(xsim[,1])), c(ourUB, rev(ourLB)),
  #        col = "skyblue", border = NA)    
  
  if (ii == 1)
    title(ylab='mRNA concentration', cex.lab = 1.5)
  
  if (ii == 4)
    title(ylab='Protein concentration', cex.lab = 1.5)
  
  if (ii == 5)
    title(xlab='Time (min)', cex.lab = 1.5)
  
  lines(xtrue[, "time"], exp(xtrue[,ii+1]),col='red', lwd=2)
  lines(xtrue[,1], ourEst, col='forestgreen', lwd=1.5)
  mtext(compnames[ii], cex=1.25)
  if (ii <= 3) points(xsim.obs$time[-1], exp(xsim.obs[-1,ii+1]), col='black', pch=16)
}

par(mar = rep(0, 4))
plot(1, type = 'n', xaxt = 'n', yaxt = 'n',
     xlab = NA, ylab = NA, frame.plot = FALSE)

legend("center", c("truth", "inferred trajectory",
                   "95% interval", "noisy observations"),
       lty = c(1, 1, 0, 0), lwd = c(2, 2, 0, 1), bty = "n",
       col = c("red", "forestgreen", NA, "black"), fill = c(0, 0, "skyblue", 0),
       border = c(0, 0, "skyblue", 0), pch = c(NA, NA, 15, 16), horiz = TRUE, cex=1.25)


## Posterior densities of parameters
par.names <- c( expression(alpha[0]), expression(alpha), "n", expression(beta))
par(mfrow=c(1,4))
for (ii in 1:4) {
  if (ii == 1) par(oma=c(0,1.5,0,0))
  par(mar = c(2.5, 2.5, 2, 0.75))
  den <- density(gpode$theta[,ii])
  plot(den, main='', xlab = '', ylab = '', type='n')
  
  value1 <- quantile(gpode$theta[,ii], 0.025)
  value2 <- quantile(gpode$theta[,ii], 0.975)
  
  l <- min(which(den$x >= value1))
  h <- max(which(den$x < value2))
  
  polygon(c(den$x[c(l, l:h, h)]),
          c(0, den$y[l:h], 0),
          col = "grey75", border=NA)
  abline(v=pram.true$theta[ii], col='red', lwd =2)
  
  lines(den)
  
  
  if (ii == 1) mtext(text='Posterior density',side=2,line=0,outer=TRUE)
  mtext(par.names[ii], cex=1.25)
}