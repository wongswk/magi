library(gpds)
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

for(kerneltype in c("compact1","rbf","matern","periodicMatern","generalMatern")){
  testpoint <- abs(rnorm(length(phitrue[[kerneltype]])+1))
  xc <- phisigllikC( testpoint, data.matrix(fn.sim[,1:2]), r, kerneltype)
  xr <- phisigllik( testpoint, fn.sim[,1:2], r, signr, TRUE, kerneltype)
  
  test_that(paste("phisigllik value -", kerneltype), {
    expect_equal(xc$value, as.numeric(xr), tolerance = 1e-4)
    expect_equal(as.numeric(xc$grad), attr(xr, "grad"), tolerance = 1e-4)
  })
  
  test_that(paste("phisigllik gradient -", kerneltype), {
    if(kerneltype=="periodicMatern"){
      skip("periodicMatern phisigllik gradient has problem for eta, needs debug")
    }
    x0 <- c(phitrue[[kerneltype]]*exp(rnorm(length(phitrue[[kerneltype]]))/10), noise)
    gradTrue <- phisigllikC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$grad
    gradNum <- c()
    for(i in 1:(length(phitrue[[kerneltype]])+1)){
      x1 = x0
      x1[i] = x1[i] + 1e-9
      gradNum[i] <- ((phisigllikC(x1, data.matrix(fn.sim[,1:2]), r, kerneltype)$value - 
                        phisigllikC(x0, data.matrix(fn.sim[,1:2]), r, kerneltype)$value)/1e-9)
    }
    expect_equal(gradNum, as.numeric(gradTrue), tolerance = 1e-4)
  })
  
  test_that(paste("phisigSample runs without error -", kerneltype), {
    if(kerneltype=="periodicMatern"){
      skip("periodicMatern phisigSample runs with error, needs debug")
    }
    phisigSample(data.matrix(fn.sim[,1:2]), r, c(phitrue[[kerneltype]], noise),
                 rep(0.03,5), 20, F, kerneltype)
  })
}

for(kerneltype in c("compact1","rbf","matern")){
  testpoint <- abs(rnorm(length(phitrue[[kerneltype]])+1))
  xc <- phisigllikC( testpoint, data.matrix(fn.sim[,1:2]), r, kerneltype)
  xc2 <- gpds:::phisigllikHard2DC( testpoint, data.matrix(fn.sim[,1:2]), r, kerneltype)
  
  test_that(paste("phisigllik returns -", kerneltype), {
    expect_equal(xc, xc2, tolerance = 1e-4)
  })
  
}
