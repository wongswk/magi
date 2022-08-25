# testhat helper functions

#' get posterior mean curve for value and derivative conditioning on
#' observed y, observed derivative dy, phi, sigma
#'
#' used to visualize Gaussian process ODE inference
#'
#' @param delta small value added to diagnal of the matrix to prevent singularity
#'
#' @noRd
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
#' @noRd
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

  infoPerRow <- 8 # FIXME grid layout not right for now
  npanel <- ceiling(ncol(infoTab)/infoPerRow)
  layout(1:npanel)
  for(i in 1:npanel){
    plot.new()
    grid::pushViewport(gridBase::baseViewports()$figure)
    gridExtra::grid.table(infoTab[,((i-1)*infoPerRow+1):min(ncol(infoTab), i*infoPerRow)])
  }
  layout(1)

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

