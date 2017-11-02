ci.list <- list()
for(it in 1:1e4){
  x <- seq(0, 20, 0.5)
  y <- rnorm(length(x), 2*x, 4)
  # plot(x, y)
  
  lambda <- 0
  beta <- list(mean = sum(x*y)/(sum(x^2) + lambda))
  sigmaSqEst <- sum((y-beta$mean*x)^2)/(length(x)-1)
  beta$sd = sqrt(sigmaSqEst/(sum(x^2) + lambda))
  
  summary(lm(y ~ 0 + x))
  
  marginalLikelihood <- function(par){
    sigmaSq <- par[1]
    lambda <- par[2]
    covKernel <- sigmaSq * (diag(length(x)) + outer(x, x)/lambda)
    -mvtnorm::dmvnorm(y, sigma = covKernel, log = TRUE)
  }
  
  optHyperparm <- optim(c(1,1), marginalLikelihood, method = "L-BFGS-B", 
                        lower = c(0,0))
  
  sigmaSqEst <- optHyperparm$par[1]
  lambda <- optHyperparm$par[2]
  
  beta <- c(mean = sum(x*y)/(sum(x^2) + lambda),
            sd = sqrt(sigmaSqEst/(sum(x^2) + lambda)))
  
  ci <- qnorm(c(0.025, 0.975), beta["mean"], beta["sd"])
  ci.list[[it]] <- ci
}

coverange <- do.call(rbind, ci.list)
mean(coverange[,1] < 2 & 2 < coverange[,2])
# ridge regression coverage is actually not bad..