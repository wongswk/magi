kerneltype <- "generalMatern"

xtime <- seq(0,10,0.01)
delta <- mean(diff(xtime))
egcov2 <- calCov(c(1,1), 
                 as.matrix(dist(xtime)),
                 -sign(outer(xtime,xtime,'-')),
                 kerneltype = kerneltype, complexity = 1)

plot(xtime, egcov2$C[,1], type="l", main=paste0(kerneltype, " C"))
plot(xtime, egcov2$Cprime[,1], type="l", main=paste0(kerneltype, " C'"))

# test with known function
plot(xtime, egcov2$Cdoubleprime[,1], type="l", main=paste0(kerneltype, " C''"))
f <- approxfun(xtime, egcov2$Cdoubleprime[,1])

intf <- approxfun(xtime, -sapply(xtime, function(ub) integrate(f, 0, ub)$value))
plot(xtime, intf(xtime), type="l")
lines(xtime, egcov2$Cprime[,1], col=2)

intintf <- approxfun(xtime, 1+sapply(xtime, function(ub) integrate(intf, 0, ub)$value))
plot(xtime, intintf(xtime), type="l")
lines(xtime, egcov2$C[,1], col=2)

# test with real application
plot(xtime, egcov2$C[,1], type="l", main=paste0(kerneltype, " C"))
f <- approxfun(xtime, egcov2$C[,1])

intf <- approxfun(xtime, -sapply(xtime, function(ub) integrate(f, 0, ub)$value))
plot(xtime, intf(xtime), type="l")

intintf <- approxfun(xtime, 1+sapply(xtime, function(ub) integrate(intf, 0, ub)$value))
plot(xtime, intintf(xtime), type="l")
lines(xtime, egcov2$C[,1], col=2)
