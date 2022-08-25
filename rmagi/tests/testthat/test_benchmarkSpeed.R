testthat::context("test benchmark speed for different likelihood")

library(testthat)
library(magi)


#### prepare for tests ####
config <- list(
  nobs = 201,
  noise = 0.1,
  seed = 123,
  bandsize = 20,
  kernel = "matern"
)

VRtrue <- read.csv(system.file("testdata/FN.csv", package="magi"))
pram.true <- list(
  abc=c(0.2,0.2,3),
  rphi=c(0.9486433, 3.2682434),
  vphi=c(1.9840824, 1.1185157),
  sigma=config$noise
)
fn.true <- VRtrue[seq(1,401,by=2),]   #### reference is 201 points
fn.true$time <- seq(0,20,0.1)
fn.sim <- fn.true

set.seed(config$seed)
fn.sim[,1:2] <- fn.sim[,1:2]+rnorm(length(unlist(fn.sim[,1:2])), sd=config$noise)
tvec.full <- fn.sim$time
fn.sim.all <- fn.sim
fn.sim[-seq(1,nrow(fn.sim), length=config$nobs),] <- NaN

foo <- outer(tvec.full, t(tvec.full),'-')[,1,]
r <- abs(foo)
r2 <- r^2
signr <- -sign(foo)

marlikmap <- list(par=c(2.314334, 1.346233, 0.622316, 2.451729, 0.084745))
cursigma <- marlikmap$par[5]

curCovV <- calCov(marlikmap$par[1:2], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)
curCovR <- calCov(marlikmap$par[3:4], r, signr, bandsize=config$bandsize, 
                  kerneltype=config$kernel)

curCovV$mu <- as.vector(fn.true[,1])  # pretend these are the means
curCovR$mu <- as.vector(fn.true[,2])

dotmu <- fnmodelODE(pram.true$abc, fn.true[,1:2]) # pretend these are the means for derivatives
curCovV$dotmu <- as.vector(dotmu[,1])  
curCovR$dotmu <- as.vector(dotmu[,2])
curCovV$tvecCovInput <- tvec.full
curCovR$tvecCovInput <- tvec.full

#### testing ####
testthat::test_that("loglik speed increase", {
  testthat::skip_on_cran()
  speedRatio <- magi:::speedbenchmarkXthetallik( data.matrix(fn.sim[,1:2]),
                                          curCovV,
                                          curCovR,
                                          cursigma,
                                          c(data.matrix(fn.true[, 1:2]), 0.2, 0.2, 3),
                                          nrep = 1e2)
  rownames(speedRatio) <- c("xthetallik_rescaled", "xthetallikBandApproxHardCode",
                            "xthetallikHardCode", "xthetallik", 
                            "xthetallik_withmu", "xthetallik_withmu2",
                            "xthetallikBandApprox", "xthetallikWithmuBand",
                            "xthetallikTwoDimension", "xthetallik in-line Band")
  testthat::expect_lt(speedRatio[2]/speedRatio[1], 0.1)
  testthat::expect_lt(speedRatio[3]/speedRatio[1], 0.5)
  # testthat::expect_lt(abs(speedRatio[3]/speedRatio[4]-1), 0.2)
  # testthat::expect_lt(abs(speedRatio[3]/speedRatio[5]-1), 0.2)
})

# lapseTime <- parallel::mclapply(1:100, function(dummy){
#   speedRatio <- speedbenchmarkXthetallik( data.matrix(fn.sim[,1:2]),
#                                           curCovV,
#                                           curCovR,
#                                           cursigma,
#                                           c(data.matrix(fn.true[, 1:2]), 0.2, 0.2, 3),
#                                           nrep = 1e3)
#   rownames(speedRatio) <- c("xthetallik_rescaled", "xthetallikBandApproxHardCode",
#                             "xthetallikHardCode", "xthetallik",
#                             "xthetallik_withmu", "xthetallik_withmu2",
#                             "xthetallikBandApprox", "xthetallikWithmuBand",
#                             "xthetallikTwoDimension", "xthetallik in-line Band")
#   speedRatio
# }, mc.cores=8)
# 
# lapseTime <- do.call(cbind, lapseTime)
# 
# rowMeans(lapseTime)/mean(lapseTime)
# t.test(lapseTime["xthetallik_withmu",] - lapseTime["xthetallik_withmu2",])
# t.test(lapseTime["xthetallik",] - lapseTime["xthetallik_withmu2",])
# t.test(lapseTime["xthetallikHardCode",] - lapseTime["xthetallik_withmu2",])

# testpoint <- abs(rnorm(5))
# microbenchmark::microbenchmark(
#   magi:::phisigllikC( testpoint, data.matrix(fn.sim[,1:2]), r, config$kernel),
#   magi:::phisigllikHard2DC( testpoint, data.matrix(fn.sim[,1:2]), r, config$kernel)
# )
