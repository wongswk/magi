library(deSolve)
library(CollocInfer)
source("R/helper/utilities.r")
source("R/helper/basic_hmc.R")
arg <- commandArgs(trailingOnly = TRUE)
lambda <- as.numeric(arg)

xsim.obs <- structure(list(time = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 
5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 
12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 
19, 19.5, 20), V = c(-0.759194784296978, 0.307652450229351, 1.76951111711703, 
1.86276223464014, 2.36260666272031, 1.94934293510063, 1.8049058730512, 
2.31411085473799, 1.27203874268536, 1.16205667831092, 0.766968584381175, 
0.942655177573709, -1.57979681438972, -1.72097792683493, -1.51357722365183, 
-1.74652934869102, -1.80665231228903, -1.36421401262715, -0.897584086467928, 
-0.20998001673955, 2.13017996012288, 1.9871890390996, 2.22364597176375, 
1.83276298302498, 1.52971311700477, 1.12786876868578, 0.908055122399705, 
0.784122302974004, 0.916015388443118, 0.452502879808095, -1.35561057200446, 
-1.9146548507507, -1.43303827258388, -1.60226388270744, -1.38937340211176, 
-1.22087533284386, -0.989167121314234, -0.00928779581803327, 
1.25484386565204, 1.75966147786693, 1.46326205875259), R = c(0.97810937587467, 
1.00832726251414, 1.08873884453796, 0.493460308463488, 0.244097745675887, 
0.192797750463905, -0.319126335971687, -0.460599665499682, -0.548733902771001, 
-0.766818946081773, -0.826547424961883, -0.76244577255139, -0.825681940592396, 
-0.411385740944421, -0.144924576130485, 0.170778189484201, 0.715204900062833, 
0.652918420188037, 0.787173269056227, 1.02645270501105, 0.844188800411126, 
0.532717423576459, 0.215860566963719, 0.122295058853959, -0.0339556924546604, 
-0.703134430649982, -0.644842746298294, -0.758511757328349, -0.957002158720152, 
-0.810190666561747, -0.894137819023785, -0.351304661067168, -0.131154633034128, 
0.365289671899918, 0.500382305671136, 0.789524449348573, 0.478327157982601, 
1.13383381845811, 0.999962167001017, 0.573197501787208, 0.164923514989322
)), row.names = c(NA, 41L), class = "data.frame")
outDir <- "~/Workspace/DynamicSys/results/ramsay/"

fhn.ode <- function(times, x, p) {
  dx <- x
  dimnames(dx) <- dimnames(x)
  dx["V"] <- p["c"] * (x["V"] - x["V"]^3 / 3 + x["R"])
  dx["R"] <- - (x["V"] - p["a"] + p["b"] * x["R"]) / p["c"]
  return(list(dx))
}
FhNvarnames <- c("V", "R")
FhNparnames <- c("a", "b", "c")

x0 <- c(-1, 1)
names(x0) <- FhNvarnames
names(xsim.obs)[-1] <- FhNvarnames

FhNpars <- c(0.2, 0.2, 3)
names(FhNpars) <- FhNparnames

fhn.fun <- function(times, x, p, more) {
  dx <- x
  dx[, "V"] <- p["c"] * (x[, "V"] - x[, "V"]^3 / 3 + x[, "R"])
  dx[, "R"] <- -(x[, "V"] - p["a"] + p["b"] * x[, "R"]) / p["c"]
  return(dx)
}

FhNtimes <- seq(0, 20, 0.1)
FhNn <- length(FhNtimes)

# set up B-splines.
FhNrange <- c(0, 20)
breaks <- seq(0, 20, 0.1)
FhNbasis <- create.bspline.basis(range = FhNrange, norder = 4, breaks = breaks)

FhNfdPar <- fdPar(FhNbasis, int2Lfd(2), 1)
if(length(lambda) == 0){
  lambda <- 0.1  # start pars at truth so this OK, Ramsey method's estimates converge as lambda --> infty  
}

# initial smooth
DEfd0 <- smooth.basis(xsim.obs$time, data.matrix(xsim.obs[, 2:3]), FhNfdPar)$fd
coefs0 <- DEfd0$coef

profile.obj <- LS.setup(pars = FhNpars, fn = fhn.fun, lambda = lambda,
                        times = FhNtimes, coefs = coefs0, basisvals = FhNbasis)

proc <- profile.obj$proc
lik <- profile.obj$lik


Ores2.2 <- Profile.LS(fhn.fun, data.matrix(xsim.obs[,-1]), xsim.obs$time, FhNpars, coefs0, FhNbasis, lambda)    # Ramsey optimization, use this as starting values for sampling

Ores2.2$pars

pdf(paste0(outDir, "illustration-ramsay-lambda-",sprintf("%03d", round(lambda)),".pdf"), width = 10, height = 9)
par(mar=c(5.1, 4.1, 4.1*1.5, 2.1)*1.2)
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=0, pch=20, xlab='time', ylab=NA)
matplot(FhNtimes, Ores2.2$coefs[c(-1,-nrow(Ores2.2$coefs)),], lty=1, add=TRUE, type="l", 
        lwd=1, col=c("thistle4", "indianred"))
matplot(xsim.obs$time, xsim.obs[,-1], type="p", col=1:(ncol(xsim.obs)-1), pch=20, xlab='time', ylab=NA, add=TRUE)
par(family="mono")
title(paste0("lambda = ", sprintf("%3d", round(lambda))), line=1, cex.main=3)
dev.off()
save.image(paste0(outDir, "illustration-ramsay-lambda-",sprintf("%03d", round(lambda)),".rda"))
