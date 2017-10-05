
loglik <- loglikRaw


# Set up for 41 points
tvec41 <- seq(0,20, by = 0.05*10)
foo <- outer(tvec41, t(tvec41),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)


# Calc log-lik of STAN output posterior means
loglik(cbind(colMeans(vdVmcurve[,seq(1,401,length=41)]),colMeans(vdRmcurve[,seq(1,401,length=41)])), colMeans(gpsmooth_ss$abc), c(colMeans(gpsmooth_ss$vphi), colMeans(gpsmooth_ss$rphi)),mean(gpsmooth_ss$sigma), fn.sim[,1:2])

# "best" log-likelihood based on truth. phi vector found by optim with other inputs set at truth
loglik( VRtrue[seq(1,401,length=41),], c(0.2,0.2,3), c(1.9840824, 1.1185157, 0.9486433, 3.2682434), 0.1,  fn.sim[,1:2])
