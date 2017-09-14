library(deSolve)

nobs.candidates <- (2:20)^2+1
noise.candidates <- seq(0.05, 1.5, 0.05)

n.iter <- 5000
burn <- n.iter * 0.3

fhn.ode <- function(times, x, p) {
  dx <- x
  dimnames(dx) <- dimnames(x)
  dx["V"] <- p["c"] * (x["V"] - x["V"]^3 / 3 + x["R"])
  dx["R"] <- - (x["V"] - p["a"] + p["b"] * x["R"]) / p["c"]
  return(list(dx))
}

FhNvarnames <- c("V", "R")
FhNparnames <- c("a", "b", "c")

resmat <- matrix(NA, nrow=570, ncol = 13)

for (nn in 1:570) {
  args <- nn
  noise <- noise.candidates[args%%length(noise.candidates)+1]
  args <- args%/%length(noise.candidates)
  nobs <- nobs.candidates[args%%length(nobs.candidates)+1]
  #print(c(noise, nobs))

  if (file.exists(paste0("../results/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda"))) {
    load(paste0("../results/C-ode-HMC-fixphi-noise",noise,"-nobs",nobs,".rda"))

    x0 <- c(-1,1)
    FhNpars <- c(0.2, 0.2, 3)
    names(x0) <- FhNvarnames
    names(FhNpars) <- FhNparnames
    solver.real.out <- lsoda(x0, times = tvec.nobs, fhn.ode, FhNpars)

    id.max <- which.max(gpode$lp__)
    x0 <- xth.formal[id.max,c(1,1+nobs)]
    FhNpars <- xth.formal[id.max,(2*nobs+1):(2*nobs+3)]
    names(x0) <- FhNvarnames
    names(FhNpars) <- FhNparnames
    solver.MAP.out <- lsoda(x0, times = tvec.nobs, fhn.ode, FhNpars)

    SSE.real <- sum(  (solver.real.out[,2:3] - fn.sim[,1:2])^2 )
    SSE.MAP <- sum(  (solver.MAP.out[,2:3] - fn.sim[,1:2])^2 )

    resmat[nn,] <- c(noise, nobs, SSE.real, SSE.MAP, xth.formal[id.max,(2*nobs+1):(2*nobs+3)],quantile(xth.formal[-burn,2*nobs+1],c(0.025,0.975)),
       quantile(xth.formal[-burn,2*nobs+2],c(0.025,0.975)), quantile(xth.formal[-burn,2*nobs+3],c(0.025,0.975)))

    print(round(resmat[nn,],3))
  }

}

write.csv(resmat, file="../results/allres.csv", row.names=F)