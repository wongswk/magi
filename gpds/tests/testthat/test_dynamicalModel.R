library(testthat)
library(gpds)

testDynamicalModel <- function(dynamicalModelList, modelName, x, theta){
  context(paste(modelName, "model, with derivatives"))
  test_that(paste(modelName, "- Dx is consistent with f"), {
    f <- dynamicalModelList$modelODE(theta, x)
    
    deltaSmall <- 1e-6
    numericalDx <- sapply(1:ncol(x), function(j){
      xnew <- x
      xnew[,j] <- xnew[,j] + deltaSmall
      (dynamicalModelList$modelODE(theta, xnew) - f)/deltaSmall
    }, simplify = "array")
    
    fdX <- dynamicalModelList$modelDx(theta, x)
    
    expect_equal(fdX, aperm(numericalDx, c(1,3,2)), tolerance = 1e-4)
  })
  
  test_that(paste(modelName, "- Dtheta is consistent with f"), {
    f <- dynamicalModelList$modelODE(theta, x)
    
    deltaSmall <- 1e-6
    numericalDtheta <- sapply(1:length(theta), function(i){
      thetaNew <- theta
      thetaNew[i] <- theta[i] + deltaSmall
      (dynamicalModelList$modelODE(thetaNew, x) - f)/deltaSmall
    }, simplify = "array")
    
    fDtheta <- dynamicalModelList$modelDtheta(theta, x)
    
    expect_equal(fDtheta, aperm(numericalDtheta, c(1,3,2)), tolerance = 1e-4)
  })
}

x <- cbind(1:4, 4:1)
theta <- c(0.2, 0.2, 3)
dynamicalModelList <- list(
  modelODE = gpds:::fnmodelODE,
  modelDtheta = gpds:::fnmodelDtheta,
  modelDx = gpds:::fnmodelDx
)
testDynamicalModel(dynamicalModelList, "FN system", x, theta)

x <- cbind(1:4, 4:1, sin(1:4))
theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
dynamicalModelList <- list(
  modelODE = gpds:::hes1modelODE,
  modelDtheta = gpds:::hes1modelDtheta,
  modelDx = gpds:::hes1modelDx
)
testDynamicalModel(dynamicalModelList, "Hes1 system", x, theta)

