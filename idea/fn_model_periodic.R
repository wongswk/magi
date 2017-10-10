x
r <- as.matrix(dist(x))
r2 <- as.matrix(dist(c(x, x+1e-9)))[1:length(x),-(1:length(x))]
sign(r2-r)

library(deSolve)
parameters <- c(a=0.2, b=0.2, c=3)
state <- c(V = -1, R = 1)

FNmodel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Vdt <- c * (V - V^3 / 3 + R)
    Rdt <- -1/c * ( V - a + b * R)
    list(c(Vdt, Rdt))
  })
}
times <- seq(0, 20, by = 0.05)
out <- ode(y = state, times = times, func = FNmodel, parms = parameters)

VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))
VRtrue$time <- seq(0,20,0.05)

matplot(VRtrue$time, VRtrue[,1:2], type='l')

matplot(out[,"time"], out[,2:3], type='l', add=TRUE, col=3:4)
head(out)
head(VRtrue)


parameters <- c(a=0.02, b=0.1, c=5)
state <- c(V = 1.2, R = 0.9)

out <- ode(y = state, times = times, func = FNmodel, parms = parameters)
matplot(out[,"time"], out[,2:3], type='l', col=3:4)
