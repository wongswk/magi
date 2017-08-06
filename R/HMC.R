numparam <- 41*2+3  # num HMC parameters
n.iter <- 5000  # number of HMC iterations
th.all <- matrix(NA,n.iter,numparam)
phisig <- matrix(NA,n.iter,5)

ref.th <- c( VRtrue[seq(1, 401, length = 41),1], VRtrue[seq(1, 401, length = 41),2], .2, .2, 3)
th.all[1,] <-  c( apply(gpfit_ss$vtrue, 2, mean), apply(gpfit_ss$rtrue, 2, mean), 1, 1, 1)
phisig[1,] <- c( apply(gpfit_ss$vphi, 2, mean), apply(gpfit_ss$rphi, 2, mean), mean(gpfit_ss$sigma))

bestCovV <- calCov( c(1.9840824, 1.1185157 ))
bestCovR <- calCov( c( 0.9486433, 3.2682434) )
loglik( VRtrue[seq(1,401,length=41),], c(0.2,0.2,3), bestCovV, bestCovR, 0.1,  fn.sim[,1:2])


loglik <- function(x, theta, CovV, CovR, sigma, y)  {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) { return(1e9)}
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  fr <- (f[,1] - CovV$mphi %*% Vsm)
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovV$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovV$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovV$C)$modulus) - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  
  
  # R
  #CovR <- calCov(phi[3:4])
  fr <- (f[,2] - CovR$mphi %*% Rsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2 - n * log(sigma), -0.5 * as.numeric(determinant(CovR$Kphi)$modulus) -0.5 * t(fr) %*% solve(CovR$Kphi) %*% fr,  -0.5 * as.numeric(determinant(CovR$C)$modulus) - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  return(ret)
  
}


xthetallik <- function(x, theta, CovV, CovR, sigma, y, grad = F)  {
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  
  if (min(theta) < 0) {
    ret <- -1e+9
    attr(ret,"grad") <- rep(1e9, 85)
    return(ret)
    
  }
  
  Vsm <- x[,1]
  Rsm <- x[,2]
  n <- nrow(x)
  
  f <- fODE(theta, x)
  res <- matrix(nrow=2,ncol=3)
  
  # V 
  #CovV <- calCov(phi[1:2])
  frV <- (f[,1] - CovV$mphi %*% Vsm)
  
  # R
  #CovR <- calCov(phi[3:4])
  frR <- (f[,2] - CovR$mphi %*% Rsm)
  
  res[1,] <- c( -0.5 * sum((Vsm - y[,1])^2) / sigma^2, -0.5 * t(frV) %*% CovV$Kinv %*% frV,   - 0.5 * t(Vsm) %*% CovV$Cinv %*% Vsm)
  res[2,] <- c( -0.5 * sum((Rsm - y[,2])^2) / sigma^2, -0.5 * t(frR) %*% CovR$Kinv %*% frR,   - 0.5 * t(Rsm) %*% CovR$Cinv %*% Rsm)
  
  ret <- sum(res)
  attr(ret,"components") <- res
  
  
  if(grad) {
    # V contrib
    Vtemp <- diag( c*(1 - x[,1]^2)) - CovV$mphi
    Rtemp <- diag( rep(c,n))
    aTemp <- rep(0,n)
    bTemp <- rep(0,n)
    cTemp <- f[,1] / c
    VC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovV$Kinv %*% frV
    
    # R contrib
    Vtemp <- diag( rep( -1/c, n) )
    Rtemp <- diag( rep( -b/c, n) ) - CovR$mphi
    aTemp <- rep(1/c,n)
    bTemp <- -Rsm/c
    cTemp <- f[,2] * (-1/c)
    RC2 <- 2 * t(cbind(Vtemp,Rtemp,aTemp,bTemp,cTemp)) %*% CovR$Kinv %*% frR
    
    C3 <- c(2 * CovV$Cinv %*% Vsm,  2 * CovR$Cinv %*% Rsm ,0,0,0)
    C1 <- c( 2 * (Vsm - y[,1]) / sigma^2 ,  2 * (Rsm - y[,2]) / sigma^2,0,0,0)
    
    attr(ret,"grad") <- c((VC2 + RC2 + C3 + C1) * (-0.5))
  }
  
  return(ret)
  
  
}

xthU <- function(q, grad=FALSE) {
  x <- cbind(q[1:41], q[42:82])
  theta <- q[83:85]
  
  #xthetallik(x,theta, c(1.9840824, 1.1185157, 0.9486433, 3.268243), 0.1, fn.sim[,1:2], grad)
  xthetallik(x,theta, curCovV, curCovR, cursigma, fn.sim[,1:2], grad)
}

curCovV <- calCov(phisig[1,1:2])
curCovR <- calCov(phisig[1,3:4])
cursigma <- phisig[1,5]
curllik <- xthU(th.all[1,])


full_llik <- c()
lliklist <- c()
lliklist[1] <- curllik
full_llik[1] <- loglik( cbind(th.all[1,1:41],th.all[1,42:82]), th.all[1,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2])
accepts <- 0
paccepts <- 0
sigmaN <- c()
deltas <- c()


for (t in 2:n.iter) {
  # Update X and theta
  foo <- basic_hmc(xthU, step=runif(1,0.0001,0.0002), nsteps= 20, initial=th.all[t-1,], return.traj = T)
  lliklist[t] <- foo$lpr
  th.all[t,] <- foo$final
  deltas[t] <- foo$delta
  accepts <- accepts + foo$acc
  
  # Update phi and sigma
  oldCovV <- curCovV
  oldCovR <- curCovR
  old_ll <- loglik( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], oldCovV, oldCovR, phisig[t-1,5], fn.sim[,1:2])
  ps_prop <- phisig[t-1,] + rnorm(5, 0, 0.01 * phisig[1,])
  propCovV <- calCov(ps_prop[1:2])
  propCovR <- calCov(ps_prop[3:4])
  prop_ll <- loglik( cbind(th.all[t,1:41], th.all[t,42:82]), th.all[t,83:85], propCovV, propCovR, ps_prop[5], fn.sim[,1:2])
  if (runif(1) < min(1,exp(prop_ll - old_ll))) {
    phisig[t,] <- ps_prop
    curCovV <- propCovV
    curCovR <- propCovR
    cursigma <- ps_prop[5]
    paccepts <- paccepts + 1
  } else {
    phisig[t,] <- phisig[t-1,]
  }

  full_llik[t] <- loglik( cbind(th.all[t,1:41],th.all[t,42:82]), th.all[t,83:85], curCovV, curCovR, cursigma,  fn.sim[,1:2])  
}  

id.best <- which.max(full_llik)
loglik( cbind(th.all[id.best,1:41],th.all[id.best,42:82]), th.all[id.best,83:85], calCov(phisig[id.best,1:2]),calCov(phisig[id.best,3:4]), phisig[id.best,5],  fn.sim[,1:2])
