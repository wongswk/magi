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


x <- cbind(1:4, 4:1, sin(1:4), cos(1:4))
theta <- c(0.014, 1.16e-9, 1.3e-9, 5e-10, 3.62e-9, 1.56e-9, 1e-8, 1e-8, 1e-8)
dynamicalModelList <- list(
  modelODE = magi:::HIVmodelODE,
  modelDtheta = magi:::HIVmodelDtheta,
  modelDx = magi:::HIVmodelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "HIV system", x, theta))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4))
theta <- c(0.022, 0.3, 0.031)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelODE,
  modelDtheta=magi:::MichaelisMentenModelDtheta,
  modelDx=magi:::MichaelisMentenModelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModel", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4))
x <- x/10
theta <- c(0.022, 0.3, 0.031)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenLogModelODE,
  modelDtheta=magi:::MichaelisMentenLogModelDtheta,
  modelDx=magi:::MichaelisMentenLogModelDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenLogModel", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4))
theta <- c(0.022, 0.3, 0.031)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelVaODE,
  modelDtheta=magi:::MichaelisMentenModelVaDtheta,
  modelDx=magi:::MichaelisMentenModelVaDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModelVa", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4), log(1:4))
theta <- c(0.022, 0.3, 0.031, 0.1, 0.5, 0.4)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelVb6pODE,
  modelDtheta=magi:::MichaelisMentenModelVb6pDtheta,
  modelDx=magi:::MichaelisMentenModelVb6pDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModelVb6p", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4), log(1:4))
theta <- c(0.022, 0.3, 0.031, 0.1)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelVb4pODE,
  modelDtheta=magi:::MichaelisMentenModelVb4pDtheta,
  modelDx=magi:::MichaelisMentenModelVb4pDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModelVb4p", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4), log(1:4))
theta <- c(0.022, 0.3, 0.031, 0.1, 0.5, 0.4, 0.01)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelVc7pODE,
  modelDtheta=magi:::MichaelisMentenModelVc7pDtheta,
  modelDx=magi:::MichaelisMentenModelVc7pDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModelVc7p", x, theta, tvec=1:4))

x <- cbind(1:4, 4:1, sin(1:4), cos(1:4), log(1:4))
theta <- c(0.022, 0.3, 0.031, 0.1, 0.2)
dynamicalModelList <- list(
  modelODE=magi:::MichaelisMentenModelVc5pODE,
  modelDtheta=magi:::MichaelisMentenModelVc5pDtheta,
  modelDx=magi:::MichaelisMentenModelVc5pDx
)
with(dynamicalModelList, testDynamicalModel(modelODE, modelDx, modelDtheta, "MichaelisMentenModelVc5p", x, theta, tvec=1:4))
