fnmodelODE <- function(theta,x,tvec) {
  V <- x[,1]
  R <- x[,2]
  
  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = theta[3] * (V - V^3 / 3.0 + R)
  result[,2] = -1.0/theta[3] * ( V - theta[1] + theta[2] * R)
  
  result
}

fnmodelDx <- function(theta,x,tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  
  V = x[,1]

  resultDx[,1,1] = theta[3] * (1 - V^2)
  resultDx[,2,1] = theta[3]
  
  resultDx[,1,2] = (-1.0 / theta[3])
  resultDx[,2,2] = ( -1.0*theta[2]/theta[3] )
  
  resultDx
}

fnmodelDtheta <- function(theta,x,tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
  V = x[,1]
  R = x[,2]
  
  resultDtheta[,3,1] = V - V^3 / 3.0 + R
  
  resultDtheta[,1,2] =  1.0 / theta[3] 
  resultDtheta[,2,2] = -R / theta[3]
  resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)
  
  resultDtheta
}


hes1modelODE  <-function(theta,x,tvec) {
   P = x[,1]
   M = x[,2]
   H = x[,3] 
  
  PMHdt <- array(0, c(nrow(x),ncol(x)))
  PMHdt[,1] = -theta[1]*P*H + theta[2]*M - theta[3]*P
  PMHdt[,2] = -theta[4]*M + theta[5]/(1+P^2)
  PMHdt[,3] = -theta[1]*P*H + theta[6]/(1+P^2) - theta[7]*H
  
  PMHdt
}


 hes1modelDx <- function(theta,x,tvec) {
   resultDx<- array(0, c(nrow(x), ncol(x), ncol(x)))
  
   P = x[,1]
   H = x[,3] 
  
  resultDx[,1,1] = -theta[1]*H - theta[3]
  resultDx[,2,1] = ( theta[2] )
  resultDx[,3,1] = -theta[1]*P
  
  resultDx[,1,2] = -2*theta[5]*P / (1.0 + P^2)^2
  resultDx[,2,2] = ( -theta[4] )
  
  resultDx[,1,3] = -theta[1]*H - 2*theta[6]*P / (1.0 + P^2)^2
  resultDx[,3,3] = -theta[1]*P - theta[7]
  
  resultDx
}


 hes1modelDtheta <- function(theta,x,tvec) {
   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
   P = x[,1]
   M = x[,2]
   H = x[,3] 
  
  resultDtheta[,1,1] = -P * H
  resultDtheta[,2,1] = M
  resultDtheta[,3,1] = -P
  
  resultDtheta[,4,2] = -M
  resultDtheta[,5,2] = 1/(1 + P^2)
  
  resultDtheta[,1,3] = -P * H
  resultDtheta[,6,3] = 1/(1 + P^2)
  resultDtheta[,7,3] = -H
  
  resultDtheta
}


 hes1logmodelODE <- function(theta,x,tvec) {
   P = exp(x[,1])
   M = exp(x[,2])
   H = exp(x[,3]) 
  
   PMHdt<- array(0, c(nrow(x),ncol(x)))
  PMHdt[,1] = -theta[1]*H + theta[2]*M/P - theta[3]
  PMHdt[,2] = -theta[4] + theta[5]/(1+P^2)/M
  PMHdt[,3] = -theta[1]*P + theta[6]/(1+P^2)/H - theta[7]
  
  PMHdt
}


 hes1logmodelDx <- function(theta,x,tvec) {
   resultDx<- array(0, c(nrow(x), ncol(x), ncol(x)))
  
   P = x[,1]
   M = x[,2] 
   H = x[,3] 
  
   expMminusP = exp(M-P)
   dP = -(1+exp(2*P))^(-2)*exp(2*P)*2
  
  resultDx[,1,1] = -theta[2]*expMminusP
  resultDx[,2,1] = theta[2]*expMminusP
  resultDx[,3,1] = -theta[1]*exp(H)
  
  resultDx[,1,2] = theta[5]*exp(-M)*dP
  resultDx[,2,2] = -theta[5]*exp(-M)/(1+exp(2*P))
  
  resultDx[,1,3] = -theta[1]*exp(P) + theta[6]*exp(-H)*dP
  resultDx[,3,3] = -theta[6]*exp(-H)/(1+exp(2*P))
  
  resultDx
}


 hes1logmodelDtheta <- function(theta,x,tvec) {
   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
   P = x[,1]
   M = x[,2]
   H = x[,3] 
  
  resultDtheta[,1,1] = -exp(H)
  resultDtheta[,2,1] = exp(M-P)
  resultDtheta[,3,1] = (-1)
  
  resultDtheta[,4,2] = (-1)
  resultDtheta[,5,2] = exp(-M)/(1+exp(2*P))
  
  resultDtheta[,1,3] = -exp(P)
  resultDtheta[,6,3] = exp(-H)/(1+exp(2*P))
  resultDtheta[,7,3] = (-1)
  
  resultDtheta
}


 HIVmodelODE <- function(theta,x,tvec) {
   T = exp(x[,1])
   Tm = exp(x[,2])
   Tw = exp(x[,3])
   Tmw = exp(x[,4])

   HIVdt<- array(0, c(nrow(x),ncol(x)))
  
  HIVdt[,1] = (theta[1] - 1e-6*theta[2]*Tm - 1e-6*theta[3]*Tw - 1e-6*theta[4]*Tmw)
  HIVdt[,2] = (theta[7] + 1e-6*theta[2]*T - 1e-6*theta[5]*Tw) + 1e-6*0.25*theta[4]*Tmw*T / Tm
  HIVdt[,3] = (theta[8] + 1e-6*theta[3]*T - 1e-6*theta[6]*Tm) + 1e-6*0.25*theta[4]*Tmw*T / Tw
  HIVdt[,4] = theta[9] + 0.5*1e-6*theta[4]*T + (1e-6*theta[5]+1e-6*theta[6])*Tw*Tm / Tmw
  
  HIVdt
}


 HIVmodelDx <- function(theta,x,tvec) {
   resultDx<- array(0, c(nrow(x), ncol(x), ncol(x)))
  
   T = exp(x[,1])
   Tm = exp(x[,2])
   Tw = exp(x[,3])
   Tmw = exp(x[,4])
  
  resultDx[,1,1] = (0)
  resultDx[,2,1] = -1e-6*theta[2]*Tm
  resultDx[,3,1] = -1e-6*theta[3]*Tw
  resultDx[,4,1] = -1e-6*theta[4]*Tmw
  
  resultDx[,1,2] = 1e-6*theta[2]*T + 1e-6*0.25*theta[4]*Tmw*T/Tm
  resultDx[,2,2] = -1e-6*0.25*theta[4]*Tmw*T / Tm
  resultDx[,3,2] = -1e-6*theta[5]*Tw
  resultDx[,4,2] = 0.25*1e-6*theta[4]*Tmw*T/Tm
  
  resultDx[,1,3] = 1e-6*theta[3]*T + 0.25*1e-6*theta[4]*Tmw*T/Tw
  resultDx[,2,3] = -1e-6*theta[6]*Tm
  resultDx[,3,3] = -1e-6*0.25*theta[4]*Tmw*T / Tw
  resultDx[,4,3] = 1e-6*0.25*theta[4]*Tmw*T/Tw
  
  resultDx[,1,4] = 1e-6*0.5*theta[4]*T
  resultDx[,2,4] = (1e-6*theta[5]+1e-6*theta[6])*Tw*Tm/Tmw
  resultDx[,3,4] = (1e-6*theta[5]+1e-6*theta[6])*Tm*Tw/Tmw
  resultDx[,4,4] = -(1e-6*theta[5]+1e-6*theta[6])*Tw*Tm/Tmw  
  
  resultDx
}


 HIVmodelDtheta <- function(theta,x,tvec) {
   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
   T = exp(x[,1])
   Tm = exp(x[,2])
   Tw = exp(x[,3])
   Tmw = exp(x[,4])
  
  
  resultDtheta[,1,1] = (1.0)
  resultDtheta[,2,1] = -1e-6*Tm 
  resultDtheta[,3,1] = -1e-6*Tw 
  resultDtheta[,4,1] = -1e-6*Tmw
  
  resultDtheta[,2,2] = 1e-6*T
  resultDtheta[,4,2] = 1e-6*0.25*Tmw*T / Tm
  resultDtheta[,5,2] = -1e-6*Tw
  resultDtheta[,7,2] = 1.0
  
  resultDtheta[,3,3] = 1e-6*T
  resultDtheta[,4,3] = 1e-6*0.25*Tmw*T / Tw
  resultDtheta[,6,3] = -1e-6*Tm
  resultDtheta[,8,3] = 1.0
  
  resultDtheta[,4,4] = 1e-6*0.5 * T
  resultDtheta[,5,4] = 1e-6*Tw * Tm / Tmw
  resultDtheta[,6,4] = 1e-6*Tw * Tm / Tmw
  resultDtheta[,9,4] = 1.0
  
  resultDtheta
}


 ptransmodelODE <- function(theta,x,tvec) {
   S = x[,1]
   dS = x[,2]
   R = x[,3]
   RS = x[,4]
   RPP = x[,5]
  
   resultdt<- array(0, c(nrow(x),ncol(x)))
  
  resultdt[,1] = -theta[1]*S - theta[2] * S * R + theta[3] * RS
  resultdt[,2] = theta[1]*S
  resultdt[,3] = -theta[2]*S*R + theta[3]*RS + theta[5] * RPP / (theta[6]+RPP)
  resultdt[,4] = theta[2]*S*R - theta[3]* RS - theta[4]*RS
  resultdt[,5] = theta[4]*RS - theta[5] * RPP / (theta[6]+RPP)
  
  resultdt
}


 ptransmodelDx <- function(theta,x,tvec) {
   resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
  
   S = x[,1]
   dS = x[,2]
   R = x[,3]
   RS = x[,4]
   RPP = x[,5]
  
  resultDx[,1,1] = -theta[1] - theta[2] * R
  resultDx[,3,1] = -theta[2] * S
  resultDx[,4,1] = (theta[3])
  
  resultDx[,1,2] = (theta[1])
  
  resultDx[,1,3] = -theta[2]*R
  resultDx[,3,3] = -theta[2]*S
  resultDx[,4,3] = (theta[3])
  resultDx[,5,3] =  theta[5] * theta[6] /  (theta[6] + RPP)^2
  
  resultDx[,1,4] = theta[2]*R
  resultDx[,3,4] = theta[2]*S
  resultDx[,4,4] = (-theta[3] - theta[4])
  
  resultDx[,4,5] = (theta[4])
  resultDx[,5,5] = -theta[5] * theta[6] /  (theta[6] + RPP)^2
  
  resultDx
}


 ptransmodelDtheta <- function(theta,x,tvec) {
   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
  
   S = x[,1]
   dS = x[,2]
   R = x[,3]
   RS = x[,4]
   RPP = x[,5]  
  
  resultDtheta[,1,1] = -S
  resultDtheta[,2,1] = -S*R
  resultDtheta[,3,1] = RS
  
  resultDtheta[,1,2] = S
  
  resultDtheta[,2,3] = -S*R
  resultDtheta[,3,3] = RS
  resultDtheta[,5,3] = RPP / (theta[6]+RPP)
  resultDtheta[,6,3] = -theta[5] * RPP / (theta[6]+RPP)^2
  
  resultDtheta[,2,4] = S*R
  resultDtheta[,3,4] = -RS
  resultDtheta[,4,4] = -RS
  
  resultDtheta[,4,5] = RS
  resultDtheta[,5,5] = - RPP / (theta[6]+RPP)
  resultDtheta[,6,5] = theta[5] * RPP / (theta[6]+RPP)^2
  
  resultDtheta
}

 
 
#' Test dynamic system model specification
#'
#' Given functions for the ODE and its gradients (with respect to the system components and parameters), verify the correctness of the gradients using numerical differentiation.
#'
#' @param modelODE function that computes the ODEs, specified with the form \eqn{f(theta, x, t)}. See examples.
#' @param modelDx function that computes the gradients of the ODEs with respect to the system components. See examples.
#' @param modelDtheta function that computes the gradients of the ODEs with respect to the parameters \eqn{\theta}. See examples.
#' @param modelName string giving a name for the model
#' @param x data matrix of system values, one column for each component, at which to test the gradients
#' @param theta vector of parameter values for \eqn{\theta}, at which to test the gradients
#' @param tvec vector of time points corresponding to the rows of \code{x}
#' 
#' @details Calls \code{\link[testthat]{test_that}} to test equality of the analytic and numeric gradients.
#' 
#' @examples 
#' # ODE system for Fitzhugh-Nagumo equations
#' fnmodelODE <- function(theta,x,t) {
#'   V <- x[,1]
#'   R <- x[,2]
#'
#'   result <- array(0, c(nrow(x),ncol(x)))
#'   result[,1] = theta[3] * (V - V^3 / 3.0 + R)
#'   result[,2] = -1.0/theta[3] * ( V - theta[1] + theta[2] * R)
#'   
#'   result
#' }
#' 
#' # Gradient with respect to system components
#' fnmodelDx <- function(theta,x,t) {
#'   resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))
#'   V = x[,1]
#'   
#'   resultDx[,1,1] = theta[3] * (1 - V^2)
#'   resultDx[,2,1] = theta[3]
#'   
#'   resultDx[,1,2] = (-1.0 / theta[3])
#'   resultDx[,2,2] = ( -1.0*theta[2]/theta[3] )
#'   
#'   resultDx
#' }
#' 
#' # Gradient with respect to parameters theta
#' fnmodelDtheta <- function(theta,x,t) {
#'   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
#'   
#'   V = x[,1]
#'   R = x[,2]
#'   
#'   resultDtheta[,3,1] = V - V^3 / 3.0 + R
#'   
#'   resultDtheta[,1,2] =  1.0 / theta[3] 
#'   resultDtheta[,2,2] = -R / theta[3]
#'   resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)
#'   
#'   resultDtheta
#' }
#' 
#' # Example incorrect gradient with respect to parameters theta
#' fnmodelDthetaWrong <- function(theta,x,t) {
#'   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
#'   
#'   V = x[,1]
#'   R = x[,2]
#'   
#'   resultDtheta[,3,1] = V - V^3 / 3.0 - R
#'   
#'   resultDtheta[,1,2] =  1.0 / theta[3] 
#'   resultDtheta[,2,2] = -R / theta[3]
#'   resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)
#'   
#'   resultDtheta
#' }
#' 
#' # Sample data for testing gradient correctness
#' tvec <- seq(0, 20, by = 0.5)
#' V <- c(-1.16, -0.18, 1.57, 1.99, 1.95, 1.85, 1.49, 1.58, 1.47, 0.96, 
#' 0.75, 0.22, -1.34, -1.72, -2.11, -1.56, -1.51, -1.29, -1.22, 
#' -0.36, 1.78, 2.36, 1.78, 1.8, 1.76, 1.4, 1.02, 1.28, 1.21, 0.04, 
#' -1.35, -2.1, -1.9, -1.49, -1.55, -1.35, -0.98, -0.34, 1.9, 1.99, 1.84)
#' R <- c(0.94, 1.22, 0.89, 0.13, 0.4, 0.04, -0.21, -0.65, -0.31, -0.65, 
#'  -0.72, -1.26, -0.56, -0.44, -0.63, 0.21, 1.07, 0.57, 0.85, 1.04, 
#'  0.92, 0.47, 0.27, 0.16, -0.41, -0.6, -0.58, -0.54, -0.59, -1.15, 
#'  -1.23, -0.37, -0.06, 0.16, 0.43, 0.73, 0.7, 1.37, 1.1, 0.85, 0.23)
#'  
#' # Correct gradients
#' testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDtheta, 
#'     "FN equations", cbind(V,R), c(.5, .6, 2), tvec)
#'     
#' # Incorrect theta gradient
#' \dontrun{
#' testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDthetaWrong, 
#'     "FN equations", cbind(V,R), c(.5, .6, 2), tvec)}
#'     
#' 
#' @export
testDynamicalModel <- function(modelODE, modelDx, modelDtheta, modelName, x, theta, tvec){
  testthat::context(paste(modelName, "model, with derivatives"))
  testthat::test_that(paste(modelName, "- Dx is consistent with f"), {
    f <- modelODE(theta, x, tvec)

    deltaSmall <- 1e-6
    numericalDx <- sapply(1:ncol(x), function(j){
      xnew <- x
      xnew[,j] <- xnew[,j] + deltaSmall
      (modelODE(theta, xnew, tvec) - f)/deltaSmall
    }, simplify = "array")

    fdX <- modelDx(theta, x, tvec)

    testthat::expect_equal(fdX, aperm(numericalDx, c(1,3,2)), tolerance = 1e-4)
  })

  testthat::test_that(paste(modelName, "- Dtheta is consistent with f"), {
    f <- modelODE(theta, x, tvec)

    deltaSmall <- 1e-6
    numericalDtheta <- sapply(1:length(theta), function(i){
      thetaNew <- theta
      thetaNew[i] <- theta[i] + deltaSmall
      (modelODE(thetaNew, x, tvec) - f)/deltaSmall
    }, simplify = "array")

    fDtheta <- modelDtheta(theta, x, tvec)

    testthat::expect_equal(fDtheta, aperm(numericalDtheta, c(1,3,2)), tolerance = 1e-4)
  })
}
