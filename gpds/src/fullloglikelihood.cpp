#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
#include "hmc.h"
#include "dynamicalSystemModels.h"

using namespace arma;


//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetaphisigmallik( const mat & xlatent, 
                       const vec & theta, 
                       const mat & phi, 
                       const vec & sigmaInput, 
                       const mat & yobs, 
                       const vec & xtimes,
                       const OdeSystem & fOdeModel) {
  mat distSigned(xtimes.size(), xtimes.size());
  for(unsigned int i = 0; i < distSigned.n_cols; i++){
    distSigned.col(i) = xtimes - xtimes(i);
  }
  
  vector<gpcov> CovAllDimensions(phi.n_cols);
  for(unsigned int j = 0; j < phi.n_cols; j++){
    CovAllDimensions[j] = generalMaternCov( phi.col(j), distSigned, 3);
  }
  
  vec xtheta = join_vert(vectorise(xlatent), theta);
  
  vec sigma(phi.n_cols);
  bool sigmaIsScaler = (sigmaInput.size() == 1);
  if (sigmaIsScaler){
    sigma.fill(as_scalar(sigmaInput));
  }else if(sigmaInput.size() == phi.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput dimension not right");
  }
  
  lp numerator = xthetallik(xtheta, CovAllDimensions, sigma, yobs, fOdeModel, false, ones(2));
  
  const mat & fderiv = fOdeModel.fOde(theta, xlatent);
  mat phiGradient(phi.n_rows, phi.n_cols);
  vec sigmaGradient(sigma.size());
  
  lp bigpost;
  for(unsigned int pDimEach = 0; pDimEach < yobs.n_cols; pDimEach++){
    const gpcov & covThisDim = CovAllDimensions[pDimEach];
    vec eigval(covThisDim.Sigma.n_cols);
    mat eigvec(covThisDim.Sigma.n_rows, covThisDim.Sigma.n_cols);
    
    eig_sym( eigval, eigvec, covThisDim.Sigma );
    vec eta = eigvec.t() * join_vert(xlatent.col(pDimEach), fderiv.col(pDimEach));
    
    bigpost.value += -0.5*covThisDim.Sigma.n_rows*log(2.0*datum::pi) - sum(log(eigval))/2.0 - 0.5*sum(square(eta) / eigval);
    vec fitLevelError = yobs.col(pDimEach) - xlatent.col(pDimEach);
    int nobs = uvec(find_finite(yobs.col(pDimEach))).size();
    fitLevelError(find_nonfinite(fitLevelError)).fill(0.0);
    bigpost.value += -0.5*log(2.0*datum::pi)*nobs - log(sigma(pDimEach))*nobs - 0.5 * sum(square(fitLevelError) / pow(sigma(pDimEach), 2));
    
    vec alpha = eigvec * (eta / eigval);
    mat facVtemp = alpha * alpha.t() - (eigvec.each_row() % (1.0 / eigval).t()) * eigvec.t();
    
    for(unsigned int i=0; i < covThisDim.dSigmadphiCube.n_slices; i++){
      phiGradient(i, pDimEach) = 0.5*accu(facVtemp % covThisDim.dSigmadphiCube.slice(i));
    }
    
    sigmaGradient(pDimEach) = sum(square(fitLevelError)) / pow(sigma(pDimEach), 3) - nobs / sigma(pDimEach);
  }
  
  if(sigmaIsScaler){
    bigpost.gradient.set_size(numerator.gradient.size() + phiGradient.size() + 1);
    bigpost.gradient(bigpost.gradient.size()-1) = sum(sigmaGradient);
  }else{
    bigpost.gradient.set_size(numerator.gradient.size() + phiGradient.size() + sigma.size());
    bigpost.gradient.subvec(numerator.gradient.size() + phiGradient.size(), bigpost.gradient.size()-1) = sigmaGradient;
  }
  bigpost.gradient.subvec(0, numerator.gradient.size()-1) = numerator.gradient;
  bigpost.gradient.subvec(numerator.gradient.size(), numerator.gradient.size() + phiGradient.size()-1) = vectorise(phiGradient);
  return bigpost;
}

//' sample from GP ODE for latent x and theta
//' @export
// [[Rcpp::export]]
Rcpp::List xthetaphisigmaSample( const arma::mat & xInitial, 
                                 const arma::vec & thetaInitial,
                                 const arma::mat & phiInitial, 
                                 const arma::vec & sigmaInitial, 
                                 const arma::mat & yobs, 
                                 const arma::vec & xtimes,
                                 const arma::vec & step,
                                 const std::string modelName = "FN",
                                 const int nsteps = 1, 
                                 const bool traj = false){
  const std::string loglikflag = "usual";
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  std::function<lp(vec)> tgt;
  if(loglikflag == "usual"){
    tgt = [&](const vec & xthetaphisigma) -> lp{
      const mat xlatent = mat(const_cast<double*>(xthetaphisigma.memptr()),
                              xInitial.n_rows, xInitial.n_cols, false, false);
      const mat theta = mat(const_cast<double*>(xthetaphisigma.memptr() + xlatent.size()),
                            thetaInitial.n_rows, thetaInitial.n_cols, false, false);
      const mat phi = mat(const_cast<double*>(xthetaphisigma.memptr() + xlatent.size() + theta.size()),
                          phiInitial.n_rows, phiInitial.n_cols, false, false);
      const mat sigma = mat(const_cast<double*>(xthetaphisigma.memptr() + xlatent.size() + theta.size() + phi.size()),
                            sigmaInitial.n_rows, sigmaInitial.n_cols, false, false);
      return xthetaphisigmallik( xlatent, theta, phi, sigma, yobs, xtimes, model);
    };
  }else if(loglikflag == "withmean"){
    throw std::runtime_error("withmean not supported yet");
  }else if(loglikflag == "band"){
    throw std::runtime_error("band not supported yet");
  }else if(loglikflag == "withmeanBand"){
    throw std::runtime_error("withmeanBand not supported yet");
  }else{
    throw std::runtime_error("loglikflag must be 'usual', 'withmean', 'band', or 'withmeanBand'");
  }
  
  unsigned int nparam = xInitial.size() + thetaInitial.size() + phiInitial.size() + sigmaInitial.size();
  vec initial(nparam);
  
  unsigned int assignmentIterator = 0;
  initial.subvec(assignmentIterator, assignmentIterator + xInitial.size() - 1) = vectorise(xInitial);
  assignmentIterator += xInitial.size();
  initial.subvec(assignmentIterator, assignmentIterator + thetaInitial.size() - 1) = vectorise(thetaInitial);
  assignmentIterator += thetaInitial.size();
  initial.subvec(assignmentIterator, assignmentIterator + phiInitial.size() - 1) = vectorise(phiInitial);
  assignmentIterator += phiInitial.size();
  initial.subvec(assignmentIterator, assignmentIterator + sigmaInitial.size() - 1) = vectorise(sigmaInitial);
  
  vec lb = ones<vec>(nparam) * (-datum::inf);
  lb.subvec(xInitial.size(), xInitial.size() + thetaInitial.size() - 1) = model.thetaLowerBound;
  lb.subvec(xInitial.size() + thetaInitial.size(), lb.size() - 1).fill(1e-2);
  
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


//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetaphisigmallikRcpp( const arma::mat & xlatent, 
                                   const arma::vec & theta, 
                                   const arma::mat & phi, 
                                   const arma::vec & sigma, 
                                   const arma::mat & yobs, 
                                   const arma::vec & xtimes,
                                   const std::string modelName = "FN"){
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  lp ret = xthetaphisigmallik(xlatent, 
                              theta, 
                              phi, 
                              sigma, 
                              yobs, 
                              xtimes,
                              model);
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}

