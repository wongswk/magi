#' get posterior mean curve for value and derivative conditioning on
#' observed y, observed derivative dy, phi, sigma
#' 
#' used to visualize Gaussian process ODE inference
#' 
#' @param delta small value added to diagnal of the matrix to prevent singularity
#' 
#' @export
getMeanDerivCurve <- function(x, y.mat, dy.mat, x.new, phi.mat, delta = 1e-9, sigma.mat, gamma.mat=NULL){
  tvec <- c(x.new,x,x.new,x)
  id.dnew <- 1:length(x.new)
  id.dobs <- 1:length(x) + length(x.new)
  id.vnew <- length(tvec)/2 + 1:length(x.new)
  id.vobs <- tail(1:length(tvec), length(x))

  foo <- outer(tvec, t(tvec),'-')[,1,]
  r <- abs(foo)
  r2 <- r^2

  signr <- -sign(foo)

  t(sapply(1:nrow(phi.mat), function(it){
    y <- y.mat[it,]
    dy <- dy.mat[it,]
    sigma <- sigma.mat[min(it, length(sigma.mat))]
    phi <- phi.mat[min(it, nrow(phi.mat)),]
    if(is.null(gamma.mat)){
      gamma <- 0
    }else{
      gamma <- gamma.mat[it]
    }

    if(is.null(sigma.mat)){
      sigma <- 0.1
    }else{
      sigma <- sigma.mat[it]
    }


    C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
    Cprime  <- signr* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
    Cdoubleprime <- (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((5*sqrt(5)*r)/(3*phi[2]^3)) - ((25*r2)/(3*phi[2]^4)))

    M <- matrix(NA, ncol=length(tvec), nrow=length(tvec))

    M[c(id.dnew,id.dobs),c(id.dnew,id.dobs)] <- Cdoubleprime[c(id.dnew,id.dobs),c(id.dnew,id.dobs)]
    M[c(id.vnew,id.vobs),c(id.vnew,id.vobs)] <- C[c(id.vnew,id.vobs),c(id.vnew,id.vobs)]

    M[c(id.dnew,id.dobs),c(id.vnew,id.vobs)] <- Cprime[c(id.dnew,id.dobs),c(id.vnew,id.vobs)]
    M[c(id.vnew,id.vobs),c(id.dnew,id.dobs)] <- t(Cprime[c(id.dnew,id.dobs),c(id.vnew,id.vobs)])

    diag(M)[id.vobs] <- diag(M)[id.vobs]+sigma^2
    diag(M)[id.dobs] <- diag(M)[id.dobs]+gamma^2

    stopifnot(is.finite(M))
    diag(M) <- diag(M)+delta

    M[c(id.vnew,id.dnew),c(id.vobs,id.dobs)]%*%solve(M[c(id.vobs,id.dobs),c(id.vobs,id.dobs)], c(y,dy))
  }))
}

#' plot posterior sample from ode inference
#' 
#' see test run cases for how to use this function
#'
#' @param filename string of pdf filename to save
#' @param fn.true VRtrue plus time and derivative
#' @param fn.sim noisy observations
#' @param gpode list of ode posterior that contains element stated in example
#' section
#' @param init list of true parameters abc and sigma
#' 
#' @importFrom gridExtra grid.table
#' @importFrom gridBase baseViewports
#' @importFrom grid pushViewport
#'
#' @export
plotPostSamples <- function(filename, fn.true, fn.sim, gpode, init, config){
  npostplot <- config$npostplot
  if(is.null(npostplot)) {
    npostplot <- 20
  }
  id.max <- c(which.max(gpode$lp__), which.max(gpode$lglik))
  if(is.null(gpode$drobs)) gpode$drobs <- t(gpode$fode[,2,])
  if(is.null(gpode$dvobs)) gpode$dvobs <- t(gpode$fode[,1,])

  infoTab <- as.data.frame(config)
  
  pdf(filename, width = 8, height = 8)
  
  if(all(c("gridExtra","gridBase") %in% rownames(installed.packages()))){
    infoPerRow <- 8 # FIXME grid layout not right for now
    npanel <- ceiling(ncol(infoTab)/infoPerRow)
    layout(1:npanel)
    for(i in 1:npanel){
      plot.new()
      grid::pushViewport(gridBase::baseViewports()$figure)
      gridExtra::grid.table(infoTab[,((i-1)*infoPerRow+1):min(ncol(infoTab), i*infoPerRow)]) 
    }
    layout(1)
  }
  
  id.plot <- seq(1,nrow(gpode$abc),length=npostplot)
  id.plot <- unique(as.integer(id.plot))
  id.plot <- unique(c(id.max, id.plot))

  vdRpostcurve <- getMeanDerivCurve(x=fn.sim$time, x.new=fn.true$time,
                                    y.mat=gpode$rtrue[id.plot,],
                                    dy.mat=gpode$drobs[id.plot,],
                                    sigma.mat = gpode$sigma[id.plot],
                                    phi.mat = gpode$rphi[id.plot,],
                                    gamma.mat=NULL)

  vdVpostcurve <- getMeanDerivCurve(x=fn.sim$time, x.new=fn.true$time,
                                    y.mat=gpode$vtrue[id.plot,],
                                    dy.mat=gpode$dvobs[id.plot,],
                                    sigma.mat = gpode$sigma[id.plot],
                                    phi.mat = gpode$vphi[id.plot,],
                                    gamma.mat=NULL)

  matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1),
          ylab="R", main="full posterior")
  points(fn.sim$time, fn.sim$Rtrue, col=2)
  matplot(fn.sim$time, t(gpode$rtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, head(t(vdRpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
  matplot(fn.sim$time, t(gpode$drobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, tail(t(vdRpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)


  matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
          ylab="V", main="full posterior")
  points(fn.sim$time, fn.sim$Vtrue, col=2)
  matplot(fn.sim$time, t(gpode$vtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, head(t(vdVpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
  matplot(fn.sim$time, t(gpode$dvobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, tail(t(vdVpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)

  layout(1:2)
  plot(gpode$lp__, type="l")
  hist(gpode$lp__)

  layout(matrix(1:4,2,byrow = TRUE))
  hist(gpode$abc[,1], main="a")
  abline(v=init$abc[1], col=2)
  abline(v=quantile(gpode$abc[,1], c(0.025, 0.975)), col=3)
  hist(gpode$abc[,2], main="b")
  abline(v=init$abc[2], col=2)
  abline(v=quantile(gpode$abc[,2], c(0.025, 0.975)), col=3)
  hist(gpode$abc[,3], main="c")
  abline(v=init$abc[3], col=2)
  abline(v=quantile(gpode$abc[,3], c(0.025, 0.975)), col=3)
  hist(gpode$sigma, main="sigma")
  abline(v=init$sigma, col=2)
  abline(v=quantile(gpode$sigma, c(0.025, 0.975)), col=3)

  plot.ts(gpode$abc[,1], main="a")
  abline(h=init$abc[1], col=2)
  plot.ts(gpode$abc[,2], main="b")
  abline(h=init$abc[2], col=2)
  plot.ts(gpode$abc[,3], main="c")
  abline(h=init$abc[3], col=2)
  plot(gpode$sigma, main="sigma",type="l")
  abline(h=init$sigma, col=2)

  hist(gpode$rphi[,1], main="phi1_R")
  abline(v=init$rphi[1], col=2)
  hist(gpode$rphi[,2], main="phi2_R")
  abline(v=init$rphi[2], col=2)
  hist(gpode$vphi[,1], main="phi1_V")
  abline(v=init$vphi[1], col=2)
  hist(gpode$vphi[,2], main="phi2_V")
  abline(v=init$vphi[2], col=2)

  plot.ts(gpode$rphi[,1], main="phi1_R")
  abline(h=init$rphi[1], col=2)
  plot.ts(gpode$rphi[,2], main="phi2_R")
  abline(h=init$rphi[2], col=2)
  plot.ts(gpode$vphi[,1], main="phi1_V")
  abline(h=init$vphi[1], col=2)
  plot.ts(gpode$vphi[,2], main="phi2_V")
  abline(h=init$vphi[2], col=2)

  layout(1)
  if(!is.null(gpode$lglik))
    plot(gpode$lglik, gpode$lp__, xlab="log likelihood from Rscript",
         ylab="log posterior from Rscript")


  matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1),
          ylab="R", main="maximum a posterior")
  points(fn.sim$time, fn.sim$Rtrue, col=2)
  matplot(fn.sim$time, t(gpode$rtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, head(t(vdRpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
  matplot(fn.sim$time, t(gpode$drobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, tail(t(vdRpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)


  matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
          ylab="V", main="maximum a posterior")
  points(fn.sim$time, fn.sim$Vtrue, col=2)
  matplot(fn.sim$time, t(gpode$vtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, head(t(vdVpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
  matplot(fn.sim$time, t(gpode$dvobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
  matplot(fn.true$time, tail(t(vdVpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)

  dev.off()
}

#' flexible plot posterior sample from ode inference
#' 
#' see test run cases for how to use this function, allow for multiple dimensions
#' of the X
#'
#' @param filename string of pdf filename to save
#' @param xtrue VRtrue plus time and derivative
#' @param xsim noisy observations
#' @param gpode list of ode posterior that contains element stated in example
#' section
#' @param param list of true parameters abc and sigma
#' 
#' @importFrom gridExtra grid.table
#' @importFrom gridBase baseViewports
#' @importFrom grid pushViewport
#'
#' @export
plotPostSamplesFlex <- function(filename, xtrue, dotxtrue, xsim, gpode, param, config, odemodel=NULL){
  npostplot <- config$npostplot
  if(is.null(npostplot)) {
    npostplot <- 20
  }
  xpostmean <- apply(gpode$xsampled, 2:3, mean)
  if(!is.null(odemodel)){
    times <- sort(unique(round(c(odemodel$times, xsim$time, xtrue[,"time"]), 7)))
    xdesolveTRUE <- deSolve::ode(y = param$x0, times = times, func = odemodel$modelODE, parms = param$theta)
    
    mapId <- which.max(gpode$lglik)
    ttheta <- gpode$theta[mapId,]
    tx0 <- gpode$xsampled[mapId,1,]
    xdesolveMAP <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    ttheta <- colMeans(gpode$theta)
    tx0 <- colMeans(gpode$xsampled[,1,])
    xdesolvePM <- deSolve::ode(y = tx0, times = times, func = odemodel$modelODE, parms = ttheta)
    
    rowId <- sapply(xsim$time, function(x) which(abs(x-times) < 1e-6))
    xdesolveTRUE.obs <- xdesolveTRUE[rowId,-1]
    xdesolveMAP.obs <- xdesolveMAP[rowId,-1]
    xdesolvePM.obs <- xdesolvePM[rowId,-1]
    
    xdesolveSamples <- parallel::mclapply(1:16, function(dummy){
      mapId <- sample(1:length(gpode$lglik), 1)
      ttheta <- gpode$theta[mapId,]
      tx0 <- gpode$xsampled[mapId,1,]
      deSolve::ode(y = tx0, times = odemodel$times, func = odemodel$modelODE, parms = ttheta)
    }, mc.cores = 8)
    
    rmseTrue <- sqrt(apply((xdesolveTRUE.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseWholeGpode <- sqrt(apply((xpostmean - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseOdeMAP <- sqrt(apply((xdesolveMAP.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    rmseOdePM <- sqrt(apply((xdesolvePM.obs - xsim[,-1])^2, 2, mean, na.rm=TRUE))
    
    rmselist <- list(
      true = paste0(round(rmseTrue, 3), collapse = "; "),
      wholeGpode = paste0(round(rmseWholeGpode, 3), collapse = "; "),
      odeMAP = paste0(round(rmseOdeMAP, 3), collapse = "; "),
      odePM = paste0(round(rmseOdePM, 3), collapse = "; ")
    )
    config$rmse <- rmselist
  }
  id.max <- which.max(gpode$lglik)
  config$noise <- paste(round(config$noise, 3), collapse = "; ")
  infoTab <- as.data.frame(config)
  
  pdf(filename, width = 8, height = 8)
  
  if(all(c("gridExtra","gridBase") %in% rownames(installed.packages()))){
    infoPerRow <- 6
    npanel <- ceiling(ncol(infoTab)/infoPerRow)
    tbls <- lapply(1:npanel, function(i){
      gridExtra::tableGrob(infoTab[,((i-1)*infoPerRow+1):min(ncol(infoTab), i*infoPerRow)]) 
    })
    do.call(gridExtra::grid.arrange, c(tbls, nrow=length(tbls)))
  }
  
  matplot(xtrue$time, xtrue[,-1], type="l", lty=1)
  matplot(xsim$time, xsim[,-1], type="p", pch=20, add=TRUE)
  
  id.plot <- seq(1,nrow(gpode$theta),length=npostplot)
  id.plot <- unique(as.integer(id.plot))
  id.plot <- unique(c(id.max, id.plot))
  
  for(j in 1:(ncol(xsim)-1)){
    matplot(xtrue$time, cbind(xtrue[,j+1], dotxtrue[,j]), type="l", lty=1, col=c(2,1),
            ylab=paste0("component-",j), main="full posterior")
    points(xsim$time, xsim[,j+1], col=2)
    matplot(xsim$time, t(gpode$xsampled[id.plot,,j]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
    matplot(xsim$time, t(gpode$xsampled[id.plot,,j]), col="skyblue",add=TRUE, type="l",lty=1)
    matplot(xsim$time, t(gpode$fode[id.plot,,j]), col="grey",add=TRUE, type="p",lty=1, pch=20)
    matplot(xsim$time, t(gpode$fode[id.plot,,j]), col="grey",add=TRUE, type="l",lty=1)
    
    if(!is.null(odemodel) && !is.null(odemodel$curCov)){
      lines(xsim$time, odemodel$curCov[[j]]$mu, col="forestgreen", lwd=2)
      lines(xsim$time, odemodel$curCov[[j]]$dotmu, col="darkgreen", lwd=2)
    }
    
    matplot(xtrue$time, cbind(xtrue[,j+1], dotxtrue[,j]), type="l", lty=1, col=c(2,1),
            ylab=paste0("component-",j), main="full posterior", add=TRUE)
  }
  
  layout(1:2)
  for(j in 1:(ncol(xsim)-1)){
    matplot(xtrue$time, cbind(xtrue[,j+1] - approx(xsim$time, odemodel$curCov[[j]]$mu, xtrue$time)$y, 
                              dotxtrue[,j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, xtrue$time)$y), 
            type="l", lty=1, col=c(2,1), ylab=paste0("component-",j), main="full posterior")
    points(xsim$time, xsim[,j+1] - odemodel$curCov[[j]]$mu, col=2)
    matplot(xsim$time, t(gpode$xsampled[id.plot,,j]) - odemodel$curCov[[j]]$mu, col="skyblue",add=TRUE, type="p",lty=1, pch=20)
    matplot(xsim$time, t(gpode$xsampled[id.plot,,j]) - odemodel$curCov[[j]]$mu, col="skyblue",add=TRUE, type="l",lty=1)
    
    matplot(xtrue$time, cbind(xtrue[,j+1] - approx(xsim$time, odemodel$curCov[[j]]$mu, xtrue$time)$y, 
                              dotxtrue[,j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, xtrue$time)$y), 
            type="l", lty=1, col=c(2,1), ylab=paste0("component-",j), main="full posterior")
    matplot(xsim$time, t(gpode$fode[id.plot,,j]) - odemodel$curCov[[j]]$dotmu, col="grey",add=TRUE, type="p",lty=1, pch=20)
    matplot(xsim$time, t(gpode$fode[id.plot,,j]) - odemodel$curCov[[j]]$dotmu, col="grey",add=TRUE, type="l",lty=1)
    
    matplot(xtrue$time, cbind(xtrue[,j+1] - approx(xsim$time, odemodel$curCov[[j]]$mu, xtrue$time)$y, 
                              dotxtrue[,j] - approx(xsim$time, odemodel$curCov[[j]]$dotmu, xtrue$time)$y), 
            type="l", lty=1, col=c(2,1), ylab=paste0("component-",j), main="full posterior", add=TRUE)
  }
  
  layout(matrix(1:4,2,byrow = TRUE))
  for(i in 1:length(param$theta)){
    hist(gpode$theta[,i], main=letters[i])
    abline(v=param$theta[i], col=2)
    abline(v=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
  }
  
  if(length(gpode$sigma) < length(gpode$lglik)){
    gpode$sigma <- matrix(gpode$sigma, nrow=1)
  }
  for(sigmaIt in 1:ncol(gpode$sigma)){
    hist(gpode$sigma[,sigmaIt], xlim=range(c(gpode$sigma[,sigmaIt], param$sigma[sigmaIt])),
         main = paste("sigma for component", sigmaIt))
    abline(v=param$sigma, col=2)
    abline(v=quantile(gpode$sigma[,sigmaIt], c(0.025, 0.975)), col=3)
  }
  
  layout(matrix(1:4,2,byrow = TRUE))
  for(i in 1:length(param$theta)){
    plot.ts(gpode$theta[,i], main=letters[i])
    abline(h=param$theta[i], col=2)
    abline(h=quantile(gpode$theta[,i], c(0.025, 0.975)), col=3)  
  }
  
  layout(1:2)
  plot(gpode$lglik, type="l")
  hist(gpode$lglik)
  
  layout(1)
  xpostmedian <- apply(gpode$xsampled, 2:3, median)
  for(j in 1:(ncol(xsim)-1)){
    matplot(xtrue$time, cbind(xtrue[,j+1], dotxtrue[,j]), type="l", lty=1, col=c(2,1),
            ylab=paste0("component-",j), main="maximum a posterior")
    points(xsim$time, xsim[,j+1], col=2)
    points(xsim$time, gpode$xsampled[id.max,,j], col="skyblue", pch=20)
    lines(xsim$time, gpode$xsampled[id.max,,j], col="skyblue")
    points(xsim$time, gpode$fode[id.max,,j], col="grey", pch=20)
    lines(xsim$time, gpode$fode[id.max,,j], col="grey")
    lines(xsim$time, xpostmean[,j], type="b", col="green")
    lines(xsim$time, xpostmedian[,j], type="b", col="brown")
    legend("topleft", c("MAP", "Post-Mean", "Post-Median"), col=c("skyblue", "green", "brown"), lty=1)
  }
  
  if(!is.null(odemodel)){
    layout(1)
    matplot(odemodel$xtrue[, "time"], odemodel$xtrue[, -1], type="l", lty=1, add=FALSE, lwd=3)
    for(xdesolveSampleEach in xdesolveSamples){
      matplot(xdesolveSampleEach[, "time"], xdesolveSampleEach[, -1], type="l", lty=4, add=TRUE)
    }
    matplot(xdesolveMAP[, "time"], xdesolveMAP[, -1], type="l", lty=3, add=TRUE, lwd=2)
    matplot(xdesolvePM[, "time"], xdesolvePM[, -1], type="l", lty=2, add=TRUE, lwd=2)
    matplot(xsim$time, xsim[,-1], type="p", col=1:(ncol(xsim)-1), pch=20, add = TRUE)
    legend("topleft", lty=c(1,3,2), legend=c("true", "MAP", "PosteriorMean"))
  }
  
  if(!is.null(odemodel) && !is.null(odemodel$curCov)){
    layout(1:(ncol(yobs)+1))
    curCov <- odemodel$curCov
    
    gpode$mxode <- sapply(1:ncol(yobs), function(j) 
      curCov[[j]]$mphi %*% (t(gpode$xsampled[,,j]) - curCov[[j]]$mu) + curCov[[j]]$dotmu, 
      simplify = "array")
    gpode$mxode <- aperm(gpode$mxode, c(2,1,3))
    gpode$odeErr <- gpode$fode - gpode$mxode
    postmeanOdeErr <- apply(gpode$odeErr, 2:3, mean)
    matplot(postmeanOdeErr, lty=1, type="l", main="postmeanOdeErr")
    for(j in 1:ncol(yobs))
      plot(xsim$time, cumsum(postmeanOdeErr[,j]), lty=1, type="l", main="cumsum postmeanOdeErr")
  }
  dev.off()
}

#' summarize and plot posterior sample from simple GP fit without ode inference
#' 
#' see test run file for example of usage
#'
#' @param filename string of pdf filename to save
#' @param fn.true VRtrue plus time
#' @param fn.sim noisy observations
#' @param gpfit list of GP posterior that contains element stated in example
#' section
#' @param init list of true parameters abc and sigma
#'
#' @export
summary.post.noODE <- function(filename, fn.true, fn.sim, gpfit, init, plotx=NULL){
  if(is.null(plotx)){
    plotx <- fn.sim$time
  }
  outIndex <- match(fn.sim$time, plotx)
  if(any(is.na(outIndex))) {
    stop("wrong plotx input: all values should be in fn.sim$time")
  }
  noeta <- is.null(gpfit$reta) || is.null(gpfit$veta)
  gpfit.post <- list()
  id.plot <- seq(1,length(gpfit$lp__),length=20)
  id.plot <- unique(as.integer(id.plot))
  id.max <- which.max(gpfit$lp__)

  init.map <- list(
    reta = if(noeta) NULL else gpfit$reta[id.max,],
    veta = if(noeta) NULL else gpfit$veta[id.max,],
    rtrue = gpfit$rtrue[id.max, outIndex],
    vtrue = gpfit$vtrue[id.max, outIndex],
    rphi = gpfit$rphi[id.max,],
    vphi = gpfit$vphi[id.max,],
    sigma = gpfit$sigma[id.max]
  )

  init.marmode <- list(
    reta = if(noeta) NULL else apply(gpfit$reta,2,mode.density),
    veta = if(noeta) NULL else apply(gpfit$veta,2,mode.density),
    rtrue = apply(gpfit$rtrue[, outIndex],2,mode.density),
    vtrue = apply(gpfit$vtrue[, outIndex],2,mode.density),
    rphi = apply(gpfit$rphi,2,mode.density),
    vphi = apply(gpfit$vphi,2,mode.density),
    sigma = mode.density(gpfit$sigma)
  )

  init.epost <- list(
    reta = if(noeta) NULL else colMeans(gpfit$reta),
    veta = if(noeta) NULL else colMeans(gpfit$veta),
    rtrue = colMeans(gpfit$rtrue[, outIndex]),
    vtrue = colMeans(gpfit$vtrue[, outIndex]),
    rphi = colMeans(gpfit$rphi),
    vphi = colMeans(gpfit$vphi),
    sigma = mean(gpfit$sigma)
  )

  pdf(filename, width = 8, height = 8)
  layout(1)
  matplot(fn.true$time, data.matrix(fn.true[,c(1:2)]), type="l", lty=1, col=c(2,1),
          ylab="R & V", main="full posterior", ylim=range(c(fn.sim$Rtrue, fn.sim$Vtrue)))
  points(fn.sim$time, fn.sim$Rtrue, col=1)
  points(fn.sim$time, fn.sim$Vtrue, col=2)
  matplot(fn.true$time, data.matrix(fn.true[,c(1:2)]), type="l", lty=1, col=c(2,1),
          ylab="R & V", main="full posterior")
  points(fn.sim$time, fn.sim$Rtrue, col=1)
  points(fn.sim$time, fn.sim$Vtrue, col=2)
  matplot(plotx, t(gpfit$rtrue[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
  matplot(plotx, t(gpfit$vtrue[id.plot,]), col="pink",add=TRUE, type="p",lty=1, pch=20)

  matplot(fn.true$time, data.matrix(fn.true[,c(1:2)]), type="l", lty=1, col=c(2,1),
          ylab="R & V", main="full posterior", add=TRUE)

  lines(plotx, colMeans(gpfit$rtrue), col=3)
  lines(plotx, colMeans(gpfit$vtrue), col=3)

  lines(plotx, gpfit$rtrue[id.max,], col=4)
  lines(plotx, gpfit$vtrue[id.max,], col=4)

  layout(matrix(1:6,3,byrow = TRUE))
  hist(gpfit$rphi[,1], probability = TRUE, breaks = 20)
  abline(v=init$rphi[1], col=2)
  abline(v=init.epost$rphi[1], col=3)
  abline(v=init.map$rphi[1], col=4)
  abline(v=init.marmode$rphi[1], col=5)
  gpfit.post[["rhpi1"]] <- plot.add.dlnorm(gpfit$rphi[,1])

  hist(gpfit$rphi[,2], probability = TRUE, breaks = 20)
  abline(v=init$rphi[2], col=2)
  abline(v=init.epost$rphi[2], col=3)
  abline(v=init.map$rphi[2], col=4)
  abline(v=init.marmode$rphi[2], col=5)
  gpfit.post[["rhpi2"]] <- plot.add.dlnorm(gpfit$rphi[,2])

  hist(gpfit$vphi[,1], probability = TRUE, breaks = 20)
  abline(v=init$vphi[1], col=2)
  abline(v=init.epost$vphi[1], col=3)
  abline(v=init.map$vphi[1], col=4)
  abline(v=init.marmode$vphi[1], col=5)
  gpfit.post[["vhpi1"]] <- plot.add.dlnorm(gpfit$vphi[,1])

  hist(gpfit$vphi[,2], probability = TRUE, breaks = 20)
  abline(v=init$vphi[2], col=2)
  abline(v=init.epost$vphi[2], col=3)
  abline(v=init.map$vphi[2], col=4)
  abline(v=init.marmode$vphi[2], col=5)
  gpfit.post[["vhpi2"]] <- plot.add.dlnorm(gpfit$vphi[,2])

  hist(gpfit$sigma, probability = TRUE, breaks = 20)
  abline(v=init$sigma, col=2)
  abline(v=init.epost$sigma, col=3)
  abline(v=init.map$sigma, col=4)
  abline(v=init.marmode$sigma, col=5)
  gpfit.post[["sigma"]] <- plot.add.dlnorm(gpfit$sigma)

  layout(matrix(1:6,3,byrow = TRUE))
  plot.ts(gpfit$rphi[,1])
  plot.ts(gpfit$rphi[,2])
  plot.ts(gpfit$vphi[,1])
  plot.ts(gpfit$vphi[,2])
  plot(gpfit$sigma,type='l')
  dev.off()

  gpfit.post <- do.call(rbind, gpfit.post)
  return(list(gpfit.post=gpfit.post, init.epost=init.epost,
              init.map=init.map, init.marmode=init.marmode))
}

#### small utility functions for visualization ####
mode.density <- function(x){
  den <- density(x)
  den$x[which.max(den$y)]
}

plot.add.dlnorm <- function(samples, col=6){
  plot.function(function(x) dlnorm(x, mean(log(samples)), sd(log(samples))),
                add=TRUE, from = min(samples), to = max(samples), n=1001, col=col)
  c(mean=mean(log(samples)), sd=sd(log(samples)))
}

#' plot posterior sample from ode inference
#' 
#' see test run cases for how to use this function
#'
#' @param filename string of pdf filename to save
#' @param fn.true VRtrue plus time and derivative
#' @param fn.sim noisy observations
#' @param gpode list of ode posterior that contains element stated in example
#' section
#' @param init list of true parameters abc and sigma
#'
#' @export
plot.post.samples <- function(filename, fn.true, fn.sim, gpode, init, npostplot = 20){
  gpconfig <- list(npostplot = 20)
  plotPostSamples(filename, fn.true, fn.sim, gpode, init, gpconfig)
  warning("plot.post.samples will be deprecated, please use plotPostSamples instead.\nTrying to avoid '.' in function names")
}
