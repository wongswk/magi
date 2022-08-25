#' The FitzHugh-Nagumo (FN) equations
#' 
#' @description
#' The classic FN equations model the spike potentials of neurons, where system components \eqn{X = (V,R)} represent the voltage and recovery variables, respectively.
#' 
#' \eqn{V} and \eqn{R} are governed by the following differential equations:
#'
#' \deqn{ \frac{dV}{dt} = c(V-\frac{V^3}{3}+R) }
#' \deqn{ \frac{dR}{dt} = -\frac{1}{c}(V-a+bR) }
#' 
#' where \eqn{\theta = (a,b,c)} are system parameters.
#' 
#' @param theta vector of parameters.
#' @param x matrix of system states (one per column) at the time points in \code{tvec}.
#' @param tvec vector of time points
#' 
#' @return
#' 
#' \code{fnmodelODE} returns an array with the values of the derivatives \eqn{\dot{X}}.
#' 
#' \code{fnmodelDx} returns a 3-D array with the values of the gradients with respect to \eqn{X}.
#' 
#' \code{fnmodelDtheta} returns a 3-D array with the values of the gradients with respect to \eqn{\theta}.
#' 
#' @examples
#' theta <- c(0.2, 0.2, 3)
#' x <- matrix(1:10, nrow = 5, ncol = 2)
#' tvec <- 1:5
#' 
#' fnmodelODE(theta, x, tvec)
#' 
#' @references
#' 
#' FitzHugh, R (1961). Impulses and Physiological States in Theoretical Models of Nerve Membrane. \emph{Biophysical Journal}, 1(6), 445–466.
#' 
#' @export
fnmodelODE <- function(theta, x, tvec) {
  V <- x[,1]
  R <- x[,2]

  result <- array(0, c(nrow(x),ncol(x)))
  result[,1] = theta[3] * (V - V^3 / 3.0 + R)
  result[,2] = -1.0/theta[3] * ( V - theta[1] + theta[2] * R)

  result
}

#' @rdname fnmodelODE
#' @export
fnmodelDx <- function(theta, x, tvec) {
  resultDx <- array(0, c(nrow(x), ncol(x), ncol(x)))

  V = x[,1]

  resultDx[,1,1] = theta[3] * (1 - V^2)
  resultDx[,2,1] = theta[3]

  resultDx[,1,2] = (-1.0 / theta[3])
  resultDx[,2,2] = ( -1.0*theta[2]/theta[3] )

  resultDx
}

#' @rdname fnmodelODE
#' @export
fnmodelDtheta <- function(theta, x, tvec) {
  resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))

  V = x[,1]
  R = x[,2]

  resultDtheta[,3,1] = V - V^3 / 3.0 + R

  resultDtheta[,1,2] =  1.0 / theta[3]
  resultDtheta[,2,2] = -R / theta[3]
  resultDtheta[,3,2] = 1.0/(theta[3]^2) * ( V - theta[1] + theta[2] * R)

  resultDtheta
}

#' Hes1 equations: oscillation of mRNA and protein levels
#' 
#' @description
#' The Hes1 equations model the oscillatory cycles of protein and messenger ribonucleic acid (mRNA) levels in cultured cells. The system components \eqn{X = (P, M, H)} represent the concentrations of protein, mRNA, and the Hes1-interacting factor that provides a negative feedback loop.
#' 
#' \eqn{P}, \eqn{M}, and \eqn{H} are governed by the following differential equations:
#'
#' \deqn{ \frac{dP}{dt} = -aPH + bM - cP }
#' \deqn{ \frac{dM}{dt} = -d_M M + \frac{e}{1 + P^2} }
#' \deqn{ \frac{dH}{dt} = -aPH + \frac{f}{1+ P^2} - gH } 
#' 
#' where \eqn{\theta = (a,b,c,d_M,e,f,g)} are system parameters.
#' 
#' @param theta vector of parameters.
#' @param x matrix of system states (one per column) at the time points in \code{tvec}.
#' @param tvec vector of time points
#' 
#' @return
#' 
#' \code{hes1modelODE} returns an array with the values of the derivatives \eqn{\dot{X}}.
#' 
#' \code{hes1modelDx} returns a 3-D array with the values of the gradients with respect to \eqn{X}.
#' 
#' \code{hes1modelDtheta} returns a 3-D array with the values of the gradients with respect to \eqn{\theta}.
#' 
#' \code{hes1logmodelODE}, \code{hes1logmodelDx}, and \code{hes1logmodelDtheta} are the log-transformed versions of \code{hes1modelODE}, \code{hes1modelDx}, and \code{hes1modelDtheta}, respectively.
#' 
#' @examples
#' theta <- c(0.022, 0.3, 0.031, 0.028, 0.5, 20, 0.3)
#' x <- matrix(1:15, nrow = 5, ncol = 3)
#' tvec <- 1:5
#' 
#' hes1modelODE(theta, x, tvec)
#' 
#' @references
#' 
#' Hirata H, Yoshiura S, Ohtsuka T, Bessho Y, Harada T, Yoshikawa K, Kageyama R (2002). Oscillatory Expression of the bHLH Factor Hes1 Regulated by a Negative Feedback Loop.
#' \emph{Science}, 298(5594), 840–843.
#' 
#' @export
hes1modelODE  <-function(theta, x, tvec) {
   P = x[,1]
   M = x[,2]
   H = x[,3]

  PMHdt <- array(0, c(nrow(x),ncol(x)))
  PMHdt[,1] = -theta[1]*P*H + theta[2]*M - theta[3]*P
  PMHdt[,2] = -theta[4]*M + theta[5]/(1+P^2)
  PMHdt[,3] = -theta[1]*P*H + theta[6]/(1+P^2) - theta[7]*H

  PMHdt
}


#' @rdname hes1modelODE
#' @export
hes1modelDx <- function(theta, x, tvec) {
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

#' @rdname hes1modelODE
#' @export
hes1modelDtheta <- function(theta, x, tvec) {
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

#' @rdname hes1modelODE
#' @export
hes1logmodelODE <- function(theta, x, tvec) {
  P = exp(x[,1])
  M = exp(x[,2])
  H = exp(x[,3])

  PMHdt<- array(0, c(nrow(x),ncol(x)))
  PMHdt[,1] = -theta[1]*H + theta[2]*M/P - theta[3]
  PMHdt[,2] = -theta[4] + theta[5]/(1+P^2)/M
  PMHdt[,3] = -theta[1]*P + theta[6]/(1+P^2)/H - theta[7]

  PMHdt
}

#' @rdname hes1modelODE
#' @export
hes1logmodelDx <- function(theta, x, tvec) {
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


#' @rdname hes1modelODE
#' @export
hes1logmodelDtheta <- function(theta, x, tvec) {
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

 
#' Protein transduction model
#' 
#' @description
#' The protein transduction equations model a biochemical reaction involving a signaling protein that degrades over time.  The system components \eqn{X = (S, S_d, R, S_R, R_{pp})} represent the levels of signaling protein, its degraded form, inactive state of \eqn{R}, \eqn{S-R} complex, and activated state of \eqn{R}.
#' 
#' \eqn{S}, \eqn{S_d}, \eqn{R}, \eqn{S_R} and \eqn{R_{pp}} are governed by the following differential equations:
#'
#' \deqn{ \frac{dS}{dt} = -k_1 \cdot S -k_2 \cdot S \cdot R + k_3 \cdot S_R }
#' \deqn{ \frac{dS_d}{dt} = k_1 \cdot S  }
#' \deqn{ \frac{dR}{dt} = -k_2 \cdot S \cdot R + k_3 \cdot S_R + \frac{V \cdot R_{pp}}{K_m + R_{pp}} }
#' \deqn{ \frac{dS_R}{dt} = k_2 \cdot S \cdot R - k_3 \cdot S_R - k_4 \cdot S_R } 
#' \deqn{ \frac{dR_{pp}}{dt} = k_4 \cdot S_R - \frac{V \cdot R_{pp}}{K_m + R_{pp}}} 
#' 
#' where \eqn{\theta = (k_1, k_2, k_3,k_4, V, K_m)} are system parameters.
#' 
#' @param theta vector of parameters.
#' @param x matrix of system states (one per column) at the time points in \code{tvec}.
#' @param tvec vector of time points
#' 
#' @return
#' 
#' \code{ptransmodelODE} returns an array with the values of the derivatives \eqn{\dot{X}}.
#' 
#' \code{ptransmodelDx} returns a 3-D array with the values of the gradients with respect to \eqn{X}.
#' 
#' \code{ptransmodelDtheta} returns a 3-D array with the values of the gradients with respect to \eqn{\theta}.
#' 
#' 
#' @examples
#' theta <- c(0.07, 0.6, 0.05, 0.3, 0.017, 0.3)
#' x <- matrix(1:25, nrow = 5, ncol = 5)
#' tvec <- 1:5
#' 
#' ptransmodelODE(theta, x, tvec)
#' 
#' @references
#' 
#' Vyshemirsky, V., & Girolami, M. A. (2008). Bayesian Ranking of Biochemical System Models. \emph{Bioinformatics}, 24(6), 833-839.
#'  
#' @export
ptransmodelODE <- function(theta, x, tvec) {
  S = x[,1]
  dS = x[,2]
  R = x[,3]
  RS = x[,4]
  RPP = x[,5]

  resultdt <- array(0, c(nrow(x),ncol(x)))

  resultdt[,1] = -theta[1]*S - theta[2] * S * R + theta[3] * RS
  resultdt[,2] = theta[1]*S
  resultdt[,3] = -theta[2]*S*R + theta[3]*RS + theta[5] * RPP / (theta[6]+RPP)
  resultdt[,4] = theta[2]*S*R - theta[3]* RS - theta[4]*RS
  resultdt[,5] = theta[4]*RS - theta[5] * RPP / (theta[6]+RPP)

  resultdt
}


#' @rdname ptransmodelODE
#' @export
ptransmodelDx <- function(theta, x, tvec) {
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


#' @rdname ptransmodelODE
#' @export
ptransmodelDtheta <- function(theta, x, tvec) {
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
#' @return A list with elements \code{testDx} and \code{testDtheta}, each with value \code{TRUE} if the corresponding gradient check passed and \code{FALSE} if not.
#'  
#' @examples 
#' # ODE system and gradients for Fitzhugh-Nagumo equations: fnmodelODE, fnmodelDx, fnmodelDtheta
#' 
#' # Example of incorrect gradient with respect to parameters theta
#' fnmodelDthetaWrong <- function(theta, x, tvec) {
#'   resultDtheta <- array(0, c(nrow(x), length(theta), ncol(x)))
#'   
#'   V = x[, 1]
#'   R = x[, 2]
#'   
#'   resultDtheta[, 3, 1] = V - V^3 / 3.0 - R
#'   
#'   resultDtheta[, 1, 2] =  1.0 / theta[3] 
#'   resultDtheta[, 2, 2] = -R / theta[3]
#'   resultDtheta[, 3, 2] = 1.0 / (theta[3]^2) * (V - theta[1] + theta[2] * R)
#'   
#'   resultDtheta
#' }
#' 
#' # Sample data for testing gradient correctness
#' data(FNdat)
#'  
#' # Correct gradients
#' testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDtheta, 
#'     "FN equations", FNdat[, c("V", "R")], c(.5, .6, 2), FNdat$time)
#'     
#' # Incorrect theta gradient (test fails)
#' testDynamicalModel(fnmodelODE, fnmodelDx, fnmodelDthetaWrong, 
#'     "FN equations", FNdat[, c("V", "R")], c(.5, .6, 2), FNdat$time)
#'     
#' 
#' @export
testDynamicalModel <- function(modelODE, modelDx, modelDtheta, modelName, x, theta, tvec){
  msg <- paste(modelName, "model, with derivatives\n")
  success <- TRUE

  f <- modelODE(theta, x, tvec)

  deltaSmall <- 1e-6
  numericalDx <- sapply(1:ncol(x), function(j){
    xnew <- x
    xnew[,j] <- xnew[,j] + deltaSmall
    (modelODE(theta, xnew, tvec) - f)/deltaSmall
  }, simplify = "array")

  fdX <- modelDx(theta, x, tvec)

  if(!all(abs(fdX - aperm(numericalDx, c(1,3,2))) < 1e-4)){
    success <- FALSE
    testDx <- FALSE
    msg <- paste0(msg, "Dx may not be consistent with f, further checking is advised\n")
  }else{
    testDx <- TRUE
  }

  #testDtheta <- testthat::test_that(paste(modelName, "- Dtheta is consistent with f"), {
  f <- modelODE(theta, x, tvec)

  deltaSmall <- 1e-6
  numericalDtheta <- sapply(1:length(theta), function(i){
    thetaNew <- theta
    thetaNew[i] <- theta[i] + deltaSmall
    (modelODE(thetaNew, x, tvec) - f)/deltaSmall
  }, simplify = "array")

  fDtheta <- modelDtheta(theta, x, tvec)

  if(!all(abs(fDtheta - aperm(numericalDtheta, c(1,3,2))) < 1e-4)){
    success <- FALSE
    msg <- paste0(msg, "Dtheta may not be consistent with f, further checking is advised\n")
    testDtheta <- FALSE
  }else{
    testDtheta <- TRUE
  }

  if(success){
    msg <- paste0(msg, "Dx and Dtheta appear to be correct\n")
  }
  cat(msg)

  list(testDx=testDx, testDtheta=testDtheta)
}
