#' Generate plots from \code{magioutput} object
#' @description Plots inferred system trajectories or diagnostic traceplots from the output of \code{MagiSolver}
#' @param x a \code{magioutput} object.
#' @param type string; the default \code{type = "traj"} plots inferred trajectories, while setting \code{type = "trace"} generates diagnostic traceplots for the MCMC samples of the parameters and log-posterior values.
#' @param obs logical; if true, points will be added on the plots for the observations when \code{type = "traj"}.
#' @param ci logical; if true, credible bands/intervals will be added to the plots.
#' @param ci.col string; color to use for credible bands.
#' @param comp.names vector of system component names, when \code{type = "traj"}. If provided, should be the same length as the number of system components in \eqn{X}.
#' @param par.names vector of parameter names, when \code{type = "trace"}. If provided, should be the same length as the number of parameters in \eqn{\theta}, or the combined length of \eqn{\theta} and \eqn{\sigma} when \code{sigma = TRUE}.
#' @param est string specifying the posterior quantity to plot as the estimate. Can be "mean", "median", "mode", or "none". Default is "mean", which plots the posterior mean of the MCMC samples.
#' @param lower the lower quantile of the credible band/interval, default is 0.025. Only used if \code{ci = TRUE}.
#' @param upper the upper quantile of the credible band/interval, default is 0.975. Only used if \code{ci = TRUE}.
#' @param sigma logical; if true, the noise levels \eqn{\sigma} will be included in the traceplots when \code{type = "trace"}.
#' @param lp logical; if true, the values of the log-posterior will be included in the traceplots when \code{type = "trace"}.
#' @param nplotcol the number of subplots per row.
#' @param ... additional arguments to \code{plot}.
#'
#' @details
#' Plots the inferred system trajectories (when \code{type = "traj"}) or diagnostic traceplots of the parameters and log-posterior (when \code{type = "trace"}) from the MCMC samples.
#' By default, the posterior mean is treated as the estimate of the trajectories and parameters (\code{est = "mean"}).
#' Alternatives are the posterior median (\code{est = "median"}, taken component-wise) and the posterior mode (\code{est = "mode"}, approximated by the MCMC sample with the highest log-posterior value).
#'
#' The default \code{type = "traj"} produces plots of the inferred trajectories and credible bands from the MCMC samples, one subplot for each system component.
#' By default, \code{lower = 0.025} and \code{upper = 0.975} produces a central 95\% credible band when \code{ci = TRUE}.
#' Adding the observed data points (\code{obs = TRUE}) can provide a visual assessment of the inferred trajectories.
#'
#' Setting \code{type = "trace"} generates diagnostic traceplots for the MCMC samples of the system parameters and the values of the log-posterior, which is a useful tool for informally assessing convergence.
#' In this case, the \code{est} and \code{ci} options add horizontal lines to the plots that indicate the estimate (in red) and credible interval (in green) for each parameter.
#'
#' @examples
#' # Set up odeModel list for the Fitzhugh-Nagumo equations
#' fnmodel <- list(
#'   fOde = fnmodelODE,
#'   fOdeDx = fnmodelDx,
#'   fOdeDtheta = fnmodelDtheta,
#'   thetaLowerBound = c(0, 0, 0),
#'   thetaUpperBound = c(Inf, Inf, Inf)
#' )
#'
#' # Example FN data
#' data(FNdat)
#' y <- setDiscretization(FNdat, by = 0.25)
#'
#' # Create magioutput from a short MagiSolver run (demo only, more iterations needed for convergence)
#' result <- MagiSolver(y, fnmodel, control = list(nstepsHmc = 20, niterHmc = 500)) 
#' 
#' # Inferred trajectories
#' plot(result, comp.names = c("V", "R"), xlab = "Time", ylab = "Level")
#' 
#' # Parameter trace plots
#' plot(result, type = "trace", par.names = c("a", "b", "c", "sigmaV", "sigmaR"), sigma = TRUE)
#' 
#' @importFrom graphics polygon
#' 
#' @export
plot.magioutput <- function(x, type = "traj", obs = TRUE, ci = TRUE, ci.col = "skyblue", comp.names, par.names, est = "mean", lower = 0.025, upper = 0.975, sigma = FALSE, lp = TRUE, nplotcol = 3, ...) {

  if (!is.magioutput(x)) 
    stop("\"x\" must be a magioutput object")

  if (est == "mode")
    lpmaxInd = which.max(x$lp)

  if (type == "traj") {
    xMean <- apply(x$xsampled, c(2, 3), mean)

    if (est == "median")
      xMed <- apply(x$xsampled, c(2, 3), median)

    if (missing(comp.names)) {
      comp.names = paste0("X[", 1:ncol(xMean), "]")
    } else if (length(comp.names) != ncol(xMean)) {
      stop(paste("vector of comp.names should be length", ncol(xMean), "to match the number of components"))
    }

    if (ci) {
      xLB <- apply(x$xsampled, c(2, 3), function(x) quantile(x, lower))
      xUB <- apply(x$xsampled, c(2, 3), function(x) quantile(x, upper))
    }

    nplotrow = ceiling(ncol(xMean) / nplotcol)
    if (ncol(xMean) < nplotcol) {
      nplotcol = ncol(xMean)
    }

    par(mfrow = c(nplotrow, nplotcol), mar = c(4, 4, 1.5, 1))
    for (i in 1:ncol(xMean)) {
      # Set up empty plot
      if (obs & ci) {
        plot(rep(x$tvec, 3), c(xLB[, i], xUB[, i], x$y[,i]), type = "n", ...)
      } else if (obs) {
        plot(rep(x$tvec, 2), c(xMean[, i], x$y[,i]), type = "n", ...)
      } else {
        plot(x$tvec, xMean[, i], type = "n", ...)
      }

      mtext(comp.names[i])
      if (ci) {
        polygon(c(x$tvec, rev(x$tvec)), c(xUB[, i], rev(xLB[, i])),
                col = ci.col, border = NA)
      }

      if (est == "mean")
        lines(x$tvec, xMean[, i], ...)

      if (est == "median")
        lines(x$tvec, xMed[, i], ...)

      if (est == "mode")
        lines(x$tvec, x$xsampled[lpmaxInd, ,i], ...)

      if (obs) {
        points(x$tvec, x$y[, i])
      }

    }
  }

  if (type == "trace") {

    if (missing(par.names)) {
      par.names = paste0("theta[", 1:ncol(x$theta), "]")
      if (sigma)
        par.names = c(par.names, paste0("sigma[", 1:ncol(x$sigma), "]"))
    } else if (length(par.names) != ncol(x$theta) + sigma * ncol(x$sigma)) {
      stop(paste("vector of par.names should be length", ncol(x$theta) + sigma * ncol(x$sigma), "to match the number of parameters"))
    }

    if (lp)
      par.names = c(par.names, "log-post")

    all_samples <- x$theta
    if (sigma)
      all_samples <- cbind(all_samples, x$sigma)
    if (lp)
      all_samples <- cbind(all_samples, x$lp)

    nplotrow = ceiling(length(par.names) / nplotcol)
    if (length(par.names) < nplotcol) {
      nplotcol = length(par.names)
    }

    par(mfrow = c(nplotrow, nplotcol), mar = c(5, 2, 1, 1))

    for (i in 1:ncol(all_samples)) {
      plot(all_samples[, i], main = par.names[i], type = "l", ylab = "", ...)

      if (est == "mean")
        abline(h = mean(all_samples[,i]), col = "red", lwd = 2)
      if (est == "median")
        abline(h = median(all_samples[,i]), col = "red", lwd = 2)
      if (est == "mode")
        abline(h = all_samples[lpmaxInd,i], col = "red", lwd = 2)

      if (ci) {
        abline(h = quantile(all_samples[,i], lower), col= "green")
        abline(h = quantile(all_samples[,i], upper), col= "green")
      }
    }


  }
}



#' flexible plot posterior sample from ode inference
#' 
#' for debugging purpose: see test run cases for how to use this function,
#' allows for multiple dimensions of X
#'
#' @param filename string of pdf filename to save
#' @param xtrue true ODE curve with time as index colume
#' @param dotxtrue true ODE derivative with time as index colume
#' @param xsim noisy observations
#' @param gpode result list of magi::MagiSolver output
#' @param param list of true parameters abc and sigma
#' @param config list of configuration settings
#' @param odemodel list of the ODE model function, times, and xtrue
#' 
#' @importFrom gridExtra grid.table
#' @importFrom gridBase baseViewports
#' @importFrom grid pushViewport
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline hist layout legend lines matplot mtext par plot.function plot.new points title
#' @importFrom stats approx density dist dlnorm dnorm fft median nobs plot.ts quantile runif sd weighted.mean
#' @importFrom utils head tail
#'
#' @noRd
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
    
    xdesolveSamples <- lapply(1:16, function(dummy){
      mapId <- sample(1:length(gpode$lglik), 1)
      ttheta <- gpode$theta[mapId,]
      tx0 <- gpode$xsampled[mapId,1,]
      deSolve::ode(y = tx0, times = odemodel$times, func = odemodel$modelODE, parms = ttheta)
    })
    
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
  
  infoPerRow <- 6
  npanel <- ceiling(ncol(infoTab)/infoPerRow)
  tbls <- lapply(1:npanel, function(i){
    gridExtra::tableGrob(infoTab[,((i-1)*infoPerRow+1):min(ncol(infoTab), i*infoPerRow)])
  })
  do.call(gridExtra::grid.arrange, c(tbls, nrow=length(tbls)))

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
  
  if(!is.null(odemodel) && !is.null(odemodel$curCov)){
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
    curCov <- odemodel$curCov
    layout(1:(length(curCov)+1))
    
    gpode$mxode <- sapply(1:length(curCov), function(j) 
      curCov[[j]]$mphi %*% (t(gpode$xsampled[,,j]) - curCov[[j]]$mu) + curCov[[j]]$dotmu, 
      simplify = "array")
    gpode$mxode <- aperm(gpode$mxode, c(2,1,3))
    gpode$odeErr <- gpode$fode - gpode$mxode
    postmeanOdeErr <- apply(gpode$odeErr, 2:3, mean)
    matplot(postmeanOdeErr, lty=1, type="l", main="postmeanOdeErr")
    for(j in 1:length(curCov))
      plot(xsim$time, cumsum(postmeanOdeErr[,j]), lty=1, type="l", main="cumsum postmeanOdeErr")
  }
  dev.off()
}
