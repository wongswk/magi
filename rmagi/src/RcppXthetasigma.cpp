#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
#include "hmc.h"
#include "band.h"
#include "dynamicalSystemModels.h"
#include "xthetasigma.h"

#include "RcppWrapper.h"

using namespace arma;


//' R wrapper for xthetasigmallik
//' Calculate the log posterior of x, theta, and sigma
//' @param xlatent the matrix of latent ODE curve
//' @param theta the parameter of ODE function
//' @param sigma the observation noise for each component of y
//' @param yobs matrix of observations
//' @param covAllDimInput list of covariance kernel objects
//' @param priorTemperatureInput the prior temperature for derivative, level, and observation, in that order
//' @param useBand boolean variable indicator to use band matrix approximation
//' @param useMean boolean variable indicator to use mean function in GP
//' @param modelName string of model name
//' export removed: function specific to pre-coded models only
//' @noRd
// [[Rcpp::export]]
Rcpp::List xthetasigmallikRcpp( const arma::mat & xlatent, 
                                const arma::vec & theta, 
                                const arma::vec & sigma, 
                                const arma::mat & yobs, 
                                const Rcpp::List & covAllDimInput,
                                const Rcpp::NumericVector & priorTemperatureInput = 1.0,
                                const bool useBand = false,
                                const bool useMean = false,
                                const std::string modelName = "FN"){
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "Hes1-log"){
    model = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  vector<gpcov> covAllDimensions(yobs.n_cols);
  for(unsigned j = 0; j < yobs.n_cols; j++){
    covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
  }
  
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  lp ret = xthetasigmallik(xlatent, 
                           theta, 
                           sigma, 
                           yobs, 
                           covAllDimensions,
                           model,
                           priorTemperature,
                           useBand,
                           useMean);
  
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}


//' MAGI posterior density
//' 
//' Computes the MAGI log-posterior value and gradient for an ODE model with the given inputs: the observations \eqn{Y}, the latent system trajectories \eqn{X},
//' the parameters \eqn{\theta}, the noise standard deviations \eqn{\sigma}, and covariance kernels.
//' 
//' @param y data matrix of observations
//' @param xlatent matrix of system trajectory values
//' @param theta vector of parameter values \eqn{\theta}
//' @param sigma vector of observation noise for each system component
//' @param covAllDimInput list of covariance kernel objects for each system component. Covariance calculations may be carried out with \code{\link{calCov}}.
//' @param odeModel list of ODE functions and inputs. See details.
//' @param priorTemperatureInput  vector of tempering factors for the GP prior, derivatives, and observations, in that order. Controls the influence of the GP prior relative to the likelihood.  Recommended values: the total number of observations divided by the total number of discretization points for the GP prior and derivatives, and 1 for the observations.
//' @param useBand logical: should the band matrix approximation be used?  If \code{TRUE}, \code{covAllDimInput} must include CinvBand, mphiBand, and KinvBand as computed by \code{\link{calCov}}.
//' 
//' @return A list with elements \code{value} for the value of the log-posterior density and \code{grad} for its gradient.
//' 
//' @examples
//' # Trajectories from the Fitzhugh-Nagumo equations
//' tvec <- seq(0, 20, 2)
//' Vtrue <- c(-1, 1.91, 1.38, -1.32, -1.5, 1.73, 1.66, 0.89, -1.82, -0.93, 1.89)
//' Rtrue <- c(1, 0.33, -0.62, -0.82, 0.5, 0.94, -0.22, -0.9, -0.08, 0.95, 0.3)
//' 
//' # Noisy observations
//' Vobs <- Vtrue + rnorm(length(tvec), sd = 0.05)
//' Robs <- Rtrue + rnorm(length(tvec), sd = 0.1)
//' 
//' # Prepare distance matrix for covariance kernel calculation
//' foo <- outer(tvec, t(tvec), '-')[, 1, ]
//' r <- abs(foo)
//' r2 <- r^2
//' signr <- -sign(foo)
//'   
//' # Choose some hyperparameter values to illustrate
//' rphi <- c(0.95, 3.27)
//' vphi <- c(1.98, 1.12)
//' phiTest <- cbind(vphi, rphi)
//' 
//' # Covariance computations
//' curCovV <- calCov(phiTest[,1], r, signr, kerneltype = "generalMatern")
//' curCovR <- calCov(phiTest[,2], r, signr, kerneltype = "generalMatern")
//' 
//' # Y and X inputs to MagiPosterior
//' yInput <- data.matrix(cbind(Vobs, Robs))
//' xlatentTest <- data.matrix(cbind(Vtrue, Rtrue))
//' 
//' # Create odeModel list for FN equations
//' fnmodel <- list(
//'   fOde = fnmodelODE,
//'   fOdeDx = fnmodelDx,
//'   fOdeDtheta = fnmodelDtheta,
//'   thetaLowerBound = c(0, 0, 0),
//'   thetaUpperBound = c(Inf, Inf, Inf)
//' )
//' 
//' MagiPosterior(yInput, xlatentTest, theta = c(0.2, 0.2, 3), sigma = c(0.05, 0.1),
//'     list(curCovV, curCovR), fnmodel)
//' 
//'
//' @export
// [[Rcpp::export]]
Rcpp::List MagiPosterior( const arma::mat & y,
                          const arma::mat & xlatent,
                          const arma::vec & theta,
                          const arma::vec & sigma,
                          const Rcpp::List & covAllDimInput,
                          const OdeSystem odeModel,
                          const Rcpp::NumericVector & priorTemperatureInput = 1.0,
                          const bool useBand = false){

    vector<gpcov> covAllDimensions(y.n_cols);
    for(unsigned j = 0; j < y.n_cols; j++){
        covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
    }

    const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);

    lp ret = xthetasigmallik(xlatent,
                             theta,
                             sigma,
                             y,
                             covAllDimensions,
                             odeModel,
                             priorTemperature,
                             useBand,
                             //useMean);
                             FALSE);

    return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                              Rcpp::Named("grad")=ret.gradient);
}

//' sample from GP ODE for latent x, theta, and sigma
//' Internal function for debugging purpose
//' @noRd
// [[Rcpp::export]]
Rcpp::List xthetasigmaSample( const arma::mat & yobs, 
                              const Rcpp::List & covList, 
                              const arma::vec & sigmaInit,
                              const arma::vec & xthetaInit, 
                              const arma::vec & step,
                              const int nsteps = 1, 
                              const bool traj = false, 
                              const std::string loglikflag = "usual",
                              const Rcpp::NumericVector & priorTemperatureInput = 1.0, 
                              const std::string modelName = "FN"){
  
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  vec sigma( yobs.n_cols);
  if(sigmaInit.size() == 1){
    sigma.fill( as_scalar( sigmaInit));
  }else if(sigmaInit.size() == yobs.n_cols){
    sigma = sigmaInit;
  }else{
    throw std::runtime_error("sigmaInit size not right");
  }
  
  vector<gpcov> covAllDimensions(covList.size());
  for(unsigned int i = 0; i < covList.size(); i++){
    covAllDimensions[i] = cov_r2cpp(covList[i]);
  }
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "Hes1-log"){
    model = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  bool useBand = false;
  if(loglikflag == "band" || loglikflag == "withmeanBand"){
    useBand = true;
  }
  bool useMean = false;
  if(loglikflag == "withmean" || loglikflag == "withmeanBand"){
    useMean = true;
  }
  
  std::function<lp(vec)> tgt = [&](const vec & xthetasigma) -> lp{
      const mat & xlatent = mat(const_cast<double*>( xthetasigma.memptr()), yobs.n_rows, yobs.n_cols, false, false);
      const vec & theta = vec(const_cast<double*>( xthetasigma.memptr() + yobs.size()), xthetaInit.size() - yobs.size(), false, false);
      const vec & sigma = vec(const_cast<double*>( xthetasigma.memptr() + yobs.size() + theta.size()), sigmaInit.size(), false, false);
      return xthetasigmallik( xlatent, 
                              theta, 
                              sigma, 
                              yobs, 
                              covAllDimensions,
                              model,
                              priorTemperature,
                              useBand,
                              useMean);
    };
  
  vec lb(xthetaInit.size() + sigmaInit.size());
  lb.subvec(0, yobs.size()-1).fill(-datum::inf);
  lb.subvec(yobs.size(), xthetaInit.size()-1) = model.thetaLowerBound;
  lb.subvec(xthetaInit.size(), xthetaInit.size() + sigmaInit.size() - 1).fill(1e-3);
  
  vec initial = join_vert(xthetaInit, sigmaInit);
  hmcstate post = basic_hmcC(tgt, initial, step, lb, {datum::inf}, nsteps, traj);
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("final")=post.final,
                                      Rcpp::Named("final.p")=post.finalp,
                                      Rcpp::Named("lpr")=post.lprvalue,
                                      Rcpp::Named("step")=post.step,
                                      Rcpp::Named("apr")=post.apr,
                                      Rcpp::Named("acc")=post.acc,
                                      Rcpp::Named("delta")=post.delta);
  if(traj){
    ret.push_back(post.trajp, "traj.p");
    ret.push_back(post.trajq, "traj.q");
    ret.push_back(post.trajH, "traj.H");
  }
  return ret;
}
