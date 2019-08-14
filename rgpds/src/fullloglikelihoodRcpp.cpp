#include "rcppgpds/classDefinition.h"
#include "rcppgpds/fullloglikelihood.h"
#include "rcppgpds/tgtdistr.h"
#include "rcppgpds/hmc.h"
#include "rcppgpds/dynamicalSystemModels.h"

using namespace arma;


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

