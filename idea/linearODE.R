library(deSolve)
library(mvtnorm)

# set up linear ODE --------------------------------------------------------

parameters <- c(a=1, b=-0.2)
state <- c(X = 1)

linearModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Xdt <- a + b * X
    list(c(Xdt))
  })
}

times <- seq(0, 20, by = 0.05)
out <- ode(y = state, times = times, func = linearModel, parms = parameters)
plot(out[,"time"], out[, "X"], type="l")


# analytical solution for ODE ---------------------------------------------

# uniform parameterization: a, b, x0
XfunAnalytic <- function(t, par){
  par <- c(-par[1]/par[2], par[3]+par[1]/par[2], par[2])
  par[1] + par[2] * exp(par[3] * t)
}

par.true <- with(as.list(c(state, parameters)), {
  c(a, b, X)
})


Xfunc <- function(t) XfunAnalytic(t, par.true)
plot.function(Xfunc, 0, 20, col=2, add=TRUE)


# simulate observations --------------------------------------------------
nobs <- 401
tI <- seq(0, 10, length=nobs)

set.seed(123)
y <- rnorm(length(tI))*0.1 + Xfunc(tI)

points(tI, y)

# frequentist inference ------------------------------------------------------

mse <- function(par){
  mean((y - XfunAnalytic(tI, par))^2)
}

frequentist <- optim(par.true, mse, method = "BFGS")
sigmaEst <- sqrt(frequentist$value)

llik <- function(par){
  sum(dnorm(y, XfunAnalytic(tI, par), sigmaEst, log = TRUE))
}

frequentist <- optim(frequentist$par, llik, method = "BFGS", hessian = TRUE, 
                     control = list(fnscale = -1))

parEst <- frequentist$par
varMat <- solve(-frequentist$hessian)

parEst
sqrt(diag(varMat))

(parEst - par.true)/sqrt(diag(varMat))

# gaussian process approach -------------------------------------------------
ndis <- 401
tdis <- seq(0, 10, length=ndis)
tAll <- sort(union(tdis, tI))

fn <- function(par) -phisigllikC( par, data.matrix(cbind(y,y)), 
                                  as.matrix(dist(tI)), "matern")$value
gr <- function(par) -as.vector(phisigllikC( par, data.matrix(cbind(y,y)), 
                                            as.matrix(dist(tI)), "matern")$grad)
marlikmap <- optim(rep(1,5), fn, gr, method="L-BFGS-B", lower = 0.0001)

phi <- marlikmap$par[1:2]
cursigma <- marlikmap$par[5]

signedDist <- outer(tAll, tAll, '-')

gpcov <- calCov(phi, abs(signedDist), -sign(signedDist))

sigmaMatPsudoInv <- matrix(0, ndis, ndis)
diag(sigmaMatPsudoInv)[match(tI, tAll)] <- 1/cursigma^2

AxInvPartConst <- gpcov$Cinv + sigmaMatPsudoInv + t(gpcov$mphi) %*% gpcov$Kinv %*% gpcov$mphi
AxInvPartConstKinvMphiSum <- gpcov$Kinv %*% gpcov$mphi
AxInvPartConstKinvMphiSum <- AxInvPartConstKinvMphiSum + t(AxInvPartConstKinvMphiSum)
KinvRowSum <- rowSums(gpcov$Kinv)
KinvSum <- sum(KinvRowSum)
muxRawPartConst <- t(gpcov$mphi) %*% KinvRowSum
yContribution <- rep(0, ndis)
yContribution[match(tI, tAll)] <- y/cursigma^2

logPosterior <- function(ab){
  a <- ab[1]
  b <- ab[2]
  
  AxInv <-  AxInvPartConst + b^2 * gpcov$Kinv - b * AxInvPartConstKinvMphiSum
  
  muxRaw <- a * muxRawPartConst - a * b * KinvRowSum + yContribution
  mux <- solve( AxInv, muxRaw)
  
  llik <- 0.5 * as.numeric(determinant(AxInv, logarithm = TRUE)$modulus)
  llik <- llik + -0.5 * (a^2 * KinvSum - sum( muxRaw * mux))
  return(as.numeric(llik))
}

logPosteriorGradient <- function(ab){
  a <- ab[1]
  b <- ab[2]
  
  AxInv <-  AxInvPartConst + b^2 * gpcov$Kinv - b * AxInvPartConstKinvMphiSum
  Ax <- solve(AxInv)
  
  bIminusM <- b*diag(ndis) - gpcov$mphi
  muTilde <- yContribution - a * t(bIminusM) %*% KinvRowSum
  
  g1 <- (sum(gpcov$Kinv %*% bIminusM %*% Ax %*%
               t(b*diag(ndis) - gpcov$mphi) %*% gpcov$Kinv) - KinvSum) * a -
    sum(t(yContribution) %*% Ax %*% t(bIminusM) %*% gpcov$Kinv)
  
  g2 <- 0.5*sum(Ax * (gpcov$Kinv%*%bIminusM + t(bIminusM)%*%gpcov$Kinv))
  g2 <- g2 - a * t(muTilde) %*% Ax %*% KinvRowSum
  g2 <- g2 - 0.5 * t(muTilde) %*% Ax %*% (gpcov$Kinv%*%bIminusM + t(bIminusM)%*%gpcov$Kinv) %*% Ax %*% muTilde
  
  return(c(g1, g2))
}

ab <- c(1, -0.2)
# testthat::expect_equal(logPosterior(ab), 12723.6919701936, tolerance = 1e-6)

map <- optim(ab, logPosterior, logPosteriorGradient, method = "BFGS",
             hessian = TRUE, control = list(fnscale = -1))
logPosteriorGradient(map$par)
mapVar <- solve(-map$hessian)
sqrt(diag(mapVar))

ab <- rnorm(2)
testthat::expect_equal(logPosteriorGradient(ab)[1],
                       (logPosterior(ab + c(1e-6,0)) - logPosterior(ab))/1e-6,
                       tolerance = 1e-3)
testthat::expect_equal(logPosteriorGradient(ab)[2],
                       (logPosterior(ab + c(0,1e-6)) - logPosterior(ab))/1e-6,
                       tolerance = 1e-3)

# compare frequentist and Gaussian process ------------------------------------

# using 2nd order approximation
parEst
sqrt(diag(varMat))

map$par
sqrt(diag(mapVar))

# using full likelihood / posterior

aCandidates <- seq(0.9,1.15, length=50)
bCandidates <- seq(-0.25, -0.15, length=50)

abGrid <- expand.grid(a=aCandidates, b=bCandidates)
llikGrid <- apply(abGrid, 1, logPosterior)
llikGrid <- llikGrid - max(llikGrid)
abGrid[,"wgts"] <- exp(llikGrid)/sum(exp(llikGrid))

histogram <- tapply(abGrid[,"wgts"], list(abGrid[,"a"], abGrid[,"b"]), identity)
filled.contour(x=aCandidates, y=bCandidates, histogram)
lattice::levelplot(wgts ~ a * b, data.frame(abGrid))

abGrid[,"freqlik"] <- dmvnorm(abGrid[,c("a","b")], parEst[1:2], varMat[1:2,1:2])
histogramFreq <- tapply(abGrid[,"freqlik"], list(abGrid[,"a"], abGrid[,"b"]), identity)

contour(x=aCandidates, y=bCandidates, histogram)
contour(x=aCandidates, y=bCandidates, histogramFreq, add = TRUE, col = "red") 

