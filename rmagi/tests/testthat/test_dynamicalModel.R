library(testthat)
library(magi)

x <- cbind(1:4, 4:1)
theta <- c(0.2, 0.2, 3)
dynamicalModelList <- list(
  modelODE = magi:::fnmodelODE,
  modelDtheta = magi:::fnmodelDtheta,
  modelDx = magi:::fnmodelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "FN system", x, theta))

x <- cbind(1:4, 4:1, sin(1:4))
theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
dynamicalModelList <- list(
  modelODE = magi:::hes1modelODE,
  modelDtheta = magi:::hes1modelDtheta,
  modelDx = magi:::hes1modelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "Hes1 system", x, theta))

x <- cbind(1:4, 4:1, sin(1:4)+2)
theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
dynamicalModelList <- list(
  modelODE = magi:::hes1logmodelODE,
  modelDtheta = magi:::hes1logmodelDtheta,
  modelDx = magi:::hes1logmodelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "Hes1-log system", x, theta))


# x <- cbind(1:4, 4:1, sin(1:4), cos(1:4))
# theta <- c(0.014, 1.16e-9, 1.3e-9, 5e-10, 3.62e-9, 1.56e-9, 1e-8, 1e-8, 1e-8)
# dynamicalModelList <- list(
#   modelODE = magi:::HIVmodelODE,
#   modelDtheta = magi:::HIVmodelDtheta,
#   modelDx = magi:::HIVmodelDx
# )
# with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "HIV system", x, theta))
