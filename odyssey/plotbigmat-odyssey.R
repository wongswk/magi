library(parallel)
library(Rcpp)

nobs.candidates <- (2:20)^2+1
noise.candidates <- seq(0.05, 1.5, 0.05)

args <- commandArgs(trailingOnly = TRUE)
if(length(args)>0){
  args <- as.numeric(args)
  noise <- noise.candidates[args%%length(noise.candidates)+1]
  args <- args%/%length(noise.candidates)
  nobs <- nobs.candidates[args%%length(nobs.candidates)+1]
  print(c(noise, nobs))
}

print(c(noise, nobs))
rdaname <- paste0("results/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda")
load(rdaname)
theta.pilot <- xth.pilot[(dim(xth.pilot)[[1]]-2):dim(xth.pilot)[[1]],,,drop=FALSE]
var.pilot <- lapply(1:nfold.pilot, function(fd) var(t(theta.pilot[,,fd])))
precisionmat.pilot <- lapply(var.pilot, solve)
consensus <- lapply(1:nfold.pilot, function(fd) precisionmat.pilot[[fd]]%*%theta.pilot[,,fd])
consensus <- solve(Reduce("+", precisionmat.pilot), Reduce("+", consensus))

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-consensus-a.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))

hist(consensus[1,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC a", probability = TRUE, 
     xlim = range(consensus[1,-(1:(ncol(consensus)/2))], gpode$abc[,1]))
hist(gpode$abc[,1], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[1], col=2)

dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-consensus-b.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))

hist(consensus[2,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC b", probability = TRUE,
     xlim = range(consensus[2,-(1:(ncol(consensus)/2))], gpode$abc[,2]))
hist(gpode$abc[,2], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[2], col=2)

dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-consensus-c.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))

hist(consensus[3,-(1:(ncol(consensus)/2))], col=rgb(0.5,0,0,0.5), 
     border=NA, main="consensus MC c", probability = TRUE,
     xlim = range(consensus[3,-(1:(ncol(consensus)/2))], gpode$abc[,3]))
hist(gpode$abc[,3], col=rgb(0,0,0.5,0.5), add=TRUE, border=NA, probability = TRUE)
abline(v=pram.true$abc[3], col=2)

dev.off()

print(c(noise, nobs))
rdaname <- paste0("results/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda")

plotname <- paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-raw.png")

load(rdaname)
init <- pram.true

id.max <- c(which.max(gpode$lp__), which.max(gpode$lglik))
if(is.null(gpode$drobs)) gpode$drobs <- t(gpode$fode[,2,])
if(is.null(gpode$dvobs)) gpode$dvobs <- t(gpode$fode[,1,])

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-raw.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
matplot(fn.true$time, data.matrix(fn.true[,c(1,2)]), type="l", lty=1:2, col=c(2,2), 
        main="observations", ylim=c(-5,5))
matplot(fn.sim$time, data.matrix(fn.sim[,c(1,2)]), type="p", col=c(2,2), pch=c(1,2),
        add=TRUE)
dev.off()

id.plot <- seq(1,nrow(gpode$abc),length=8)
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

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-postR.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))

matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1), 
        ylab="R", main="full posterior R")
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(gpode$rtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdRpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpode$drobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdRpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-postV.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
        ylab="V", main="full posterior V")
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpode$vtrue[id.plot,]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdVpostcurve),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpode$dvobs[id.plot,]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdVpostcurve),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-hista.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
hist(gpode$abc[,1], main="a")
abline(v=init$abc[1], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-histb.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
hist(gpode$abc[,2], main="b")
abline(v=init$abc[2], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-histc.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
hist(gpode$abc[,3], main="c")
abline(v=init$abc[3], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-tsa.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
plot.ts(gpode$abc[,1], main="a")
abline(h=init$abc[1], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-tsb.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
plot.ts(gpode$abc[,2], main="b")
abline(h=init$abc[2], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-tsc.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
plot.ts(gpode$abc[,3], main="c")
abline(h=init$abc[3], col=2)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-mapR.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
matplot(fn.true$time, data.matrix(fn.true[,c(2,5)]), type="l", lty=1, col=c(2,1),
        ylab="R", main="maximum a posterior R")
points(fn.sim$time, fn.sim$Rtrue, col=2)
matplot(fn.sim$time, t(gpode$rtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdRpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpode$drobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdRpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)
dev.off()

png(paste0("plot/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,"-mapV.png"), 
    width = 3, height = 3, units = "in", res=40)
par(mar=c(2,1,1,1))
matplot(fn.true$time, data.matrix(fn.true[,c(1,4)]), type="l", lty=1, col=c(2,1),
        ylab="V", main="maximum a posterior V")
points(fn.sim$time, fn.sim$Vtrue, col=2)
matplot(fn.sim$time, t(gpode$vtrue[id.max,,drop=FALSE]), col="skyblue",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, head(t(vdVpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="skyblue",add=TRUE, type="l",lty=1)
matplot(fn.sim$time, t(gpode$dvobs[id.max,,drop=FALSE]), col="grey",add=TRUE, type="p",lty=1, pch=20)
matplot(fn.true$time, tail(t(vdVpostcurve[1:length(id.max),,drop=FALSE]),nrow(fn.true)), col="grey",add=TRUE, type="l",lty=1)

dev.off()

save.image(rdaname)