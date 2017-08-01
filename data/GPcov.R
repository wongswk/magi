tvec <- seq(0,20, by = 0.05)
foo <- outer(tvec, t(tvec),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

## Matern covariance has 2 params

phi <- c(1,1)

C <- phi[1] * (1 + ((sqrt(5)*r)/phi[2]) + ((5*r2)/(3*phi[2]^2))) * exp((-sqrt(5)*r)/phi[2])
Cprime  <- signr* (phi[1] * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3)))
Cdoubleprime <- (-phi[1] * (sqrt(5)/phi[2]) * exp((-sqrt(5)*r)/phi[2])) * (((5*r)/(3*phi[2]^2)) + ((5*sqrt(5)*r2)/(3*phi[2]^3))) + (phi[1]*exp((-sqrt(5)*r)/phi[2])) * ((5/(3*phi[2]^2)) + ((10*sqrt(5)*r)/(3*phi[2]^3)))

Cinv <- solve(C)
mphi <-  Cprime*Cinv
Kphi <- Cdoubleprime - (Cprime %*% Cinv %*% t(Cprime))
