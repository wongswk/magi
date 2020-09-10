context("phisigloocv")
VRtrue <- read.csv(system.file("testdata/FN.csv", package="gpds"))

phitrue <- list(
  compact1 = c(2.618, 6.381, 0.152, 9.636),
  rbf = c(0.838, 0.307, 0.202, 0.653),
  matern = c(2.04, 1.313, 0.793, 3.101),
  periodicMatern = c(2.04, 1.313, 9, 0.793, 3.101, 9),
  generalMatern = c(2.04, 1.313, 0.793, 3.101)
)

nobs <- 41
noise <- 0.5

fn.true <- VRtrue
fn.true$time <- seq(0,20,0.05)
fn.sim <- fn.true
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=noise)
fn.sim <- fn.sim[seq(1,nrow(fn.sim), length=nobs),]

tvec.nobs <- fn.sim$time
foo <- outer(tvec.nobs, t(tvec.nobs),'-')[,1,]
r <- abs(foo)
signr <- -sign(foo)

for(kerneltype in c("compact1","rbf","matern","generalMatern")){
  testpoint <- abs(rnorm(length(phitrue[[kerneltype]])+1))
  xc <- phisigloocvllikC( testpoint, data.matrix(fn.sim[,1:2]), r, kerneltype)
  
  test_that(paste("phisigloocvllikC gradient -", kerneltype), {
    x0 <- c(phitrue[[kerneltype]]*exp(rnorm(length(phitrue[[kerneltype]]))/10), noise)
    gradTrue <- phisigloocvllikC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$grad
    gradNum <- c()
    for(i in 1:(length(phitrue[[kerneltype]])+1)){
      x1 = x0
      x1[i] = x1[i] + 1e-9
      gradNum[i] <- ((phisigloocvllikC(x1, data.matrix(fn.sim[,1:2]), r, kerneltype)$value - 
                        phisigloocvllikC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$value)/1e-9)
    }
    expect_equal(gradNum, as.numeric(gradTrue), tolerance = 1e-3)
  })
  
}

for(kerneltype in c("compact1","rbf","matern","generalMatern")){
  testpoint <- abs(rnorm(length(phitrue[[kerneltype]])+1))
  xc <- phisigloocvmseC( testpoint, data.matrix(fn.sim[,1:2]), r, kerneltype)
  testpoint2 <- testpoint
  multiplyer <- rexp(1)
  testpoint[c(1,3)] <- multiplyer*testpoint[c(1,3)]
  testpoint[5] <- sqrt(multiplyer)*testpoint[5]
  xc2 <- phisigloocvmseC( testpoint2, data.matrix(fn.sim[,1:2]), r, kerneltype)
  test_that("phisigloocvmseC doesn't know variance", expect_equal(xc, xc2))
  
  test_that(paste("phisigloocvmseC gradient -", kerneltype), {
    x0 <- c(phitrue[[kerneltype]]*exp(rnorm(length(phitrue[[kerneltype]]))/10), noise)
    gradTrue <- phisigloocvmseC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$grad
    gradNum <- c()
    for(i in 1:(length(phitrue[[kerneltype]])+1)){
      x1 = x0
      x1[i] = x1[i] + 1e-9
      gradNum[i] <- ((phisigloocvmseC(x1, data.matrix(fn.sim[,1:2]), r, kerneltype)$value - 
                        phisigloocvmseC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$value)/1e-9)
    }
    expect_equal(gradNum, as.numeric(gradTrue), tolerance = 1e-3)
  })
  
}
