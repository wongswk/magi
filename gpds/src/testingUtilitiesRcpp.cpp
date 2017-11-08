#include <RcppArmadillo.h>
#include <chrono>

#include "hmc.h"
#include "classDefinition.h"
#include "wrapper.h"
#include "dynamicalSystemModels.h"
#include "band.h"

using namespace arma;
// using namespace Rcpp;

//' R wrapper for basic_hmcC for normal distribution
//' 
// [[Rcpp::export]]
Rcpp::List hmcNormal(arma::vec initial, arma::vec step, arma::vec lb, arma::vec ub,
               int nsteps = 1, bool traj = false){
  // cout << "step.size() = " << step.size() << "\tinitial.size() = " << initial.size() << endl;
  // cout << "step.n_rows = " << step.n_rows << "\tinitial.n_rows = " << initial.n_rows << endl;
  hmcstate post = basic_hmcC(lpnormal, initial, step, 
                             lb, 
                             ub, 
                             nsteps, traj);
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

// old log likelihood for latent states and ODE theta conditional on phi sigma
// 
// use for benchmarking about speed
lp xthetallikHardCode( const vec & xtheta, 
                       const gpcov & CovV, 
                       const gpcov & CovR, 
                       const double & sigma, 
                       const mat & yobs, 
                       const std::function<mat (vec, mat)> & fODE) {
  int n = (xtheta.size() - 3)/2;
  const vec & theta = xtheta.subvec(xtheta.size() - 3, xtheta.size() - 1);
  lp ret;
  
  if (min(theta) < 0) {
    ret.value = -1e+9;
    ret.gradient = zeros<vec>(2*n);
    ret.gradient.subvec(xtheta.size() - 3, xtheta.size() - 1).fill(1e9);
    return ret;
  }
  
  const vec & Vsm = xtheta.subvec(0, n - 1);
  const vec & Rsm = xtheta.subvec(n, 2*n - 1);
  
  const mat & fderiv = fODE(theta, join_horiz(Vsm, Rsm));
  mat res(2,3);
  
  // V 
  vec frV = (fderiv.col(0) - CovV.mphi * Vsm); // n^2 operation
  vec VsmCTrans = CovV.CeigenVec.t() * Vsm;
  // vec frV = fderiv.col(0) - CovV.mphiLeftHalf * (VsmCTrans % CovV.Ceigen1over);
  vec frVKTrans = CovV.KeigenVec.t() * frV;
  vec fitLevelErrorV = Vsm - yobs.col(0);
  fitLevelErrorV(find_nonfinite(fitLevelErrorV)).fill(0.0);
  res(0,0) = -0.5 * sum(square( fitLevelErrorV )) / pow(sigma,2);
  res(0,1) = -0.5 * sum( square(frVKTrans) % CovV.Keigen1over);
  res(0,2) = -0.5 * sum( square(VsmCTrans) % CovV.Ceigen1over);
  
  // R
  vec frR = (fderiv.col(1) - CovR.mphi * Rsm); // n^2 operation
  vec RsmCTrans = CovR.CeigenVec.t() * Rsm;
  // vec frR = fderiv.col(1) - CovR.mphiLeftHalf * (RsmCTrans % CovR.Ceigen1over);
  vec frRKTrans = CovR.KeigenVec.t() * frR;
  vec fitLevelErrorR = Rsm - yobs.col(1);
  fitLevelErrorR(find_nonfinite(fitLevelErrorR)).fill(0.0);
  
  res(1,0) = -0.5 * sum(square( fitLevelErrorR )) / pow(sigma,2);
  res(1,1) = -0.5 * sum( square(frRKTrans) % CovR.Keigen1over);
  res(1,2) = -0.5 * sum( square(RsmCTrans) % CovR.Ceigen1over);
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = -CovV.mphi;
  Vtemp.diag() += theta(2)*(1 - square(Vsm));
  
  vec KinvFrV = (CovV.KeigenVec * (frVKTrans % CovV.Keigen1over));
  vec abcTemp = zeros<vec>(3);
  abcTemp(2) = sum(KinvFrV % fderiv.col(0)) / theta(2);
  vec VC2 =  2.0 * join_vert(join_vert( Vtemp.t()*KinvFrV, // n^2 operation
                                        theta(2) * KinvFrV ),
                                        abcTemp );
  
  
  // R contrib
  mat Rtemp = -CovR.mphi;
  Rtemp.diag() -= theta(1)/theta(2);
  
  vec KinvFrR = (CovR.KeigenVec * (frRKTrans % CovR.Keigen1over));
  abcTemp.fill(0);
  abcTemp(0) = sum(KinvFrR) / theta(2);
  abcTemp(1) = -sum(Rsm % KinvFrR) / theta(2);
  abcTemp(2) = -sum(fderiv.col(1) % KinvFrR) / theta(2);
  vec RC2 = 2.0 * join_vert(join_vert( -KinvFrR / theta(2),
                                       Rtemp.t() * KinvFrR), // n^2 operation
                                       abcTemp );
  
  vec C3 = join_vert(join_vert( 2.0 * CovV.CeigenVec * (VsmCTrans % CovV.Ceigen1over),  
                                2.0 * CovR.CeigenVec * (RsmCTrans % CovR.Ceigen1over) ), 
                                zeros<vec>(theta.size()));
  vec C1 = join_vert(join_vert( 2.0 * fitLevelErrorV / pow(sigma,2) ,  
                                2.0 * fitLevelErrorR / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}

// log likelihood for latent states and ODE theta conditional on phi sigma
// 
// use for testing that zero mean can be used with a shifted ODE
lp xthetallik_withmu2( const vec & xtheta, 
                       const gpcov & CovV, 
                       const gpcov & CovR, 
                       const double & sigma, 
                       const mat & yobs, 
                       const OdeSystem & fOdeModel) {
  int n = (xtheta.size() - 3)/2;
  vec xthetaShifted = xtheta;
  xthetaShifted.subvec(0, n - 1) -= CovV.mu;
  xthetaShifted.subvec(n, 2*n - 1) -= CovR.mu;
  
  mat yobsShifted = yobs;
  yobsShifted.col(0) -= CovV.mu;
  yobsShifted.col(1) -= CovR.mu;
  
  OdeSystem fOdeModelShifted;
  
  fOdeModelShifted.fOde = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> mat{
      return fOdeModel.fOde(theta, x+join_horiz(CovV.mu, CovR.mu)) - join_horiz(CovV.dotmu, CovR.dotmu);
  };
  
  fOdeModelShifted.fOdeDx = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> cube{ 
    return fOdeModel.fOdeDx(theta, x+join_horiz(CovV.mu, CovR.mu));
  };
  
  fOdeModelShifted.fOdeDtheta = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> cube{ 
    return fOdeModel.fOdeDtheta(theta, x+join_horiz(CovV.mu, CovR.mu));
  };
  
  return xthetallik(xthetaShifted, CovV, CovR, sigma, yobsShifted, fOdeModelShifted);
}

//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallik_withmu2C(const arma::mat & yobs, 
                              const Rcpp::List & covVr, 
                              const Rcpp::List & covRr, 
                              const double & sigma, 
                              const arma::vec & initial){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta);
  lp ret = xthetallik_withmu2(initial, covV, covR, sigma, yobs, fnmodel);
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}


lp xthetallikBandApproxHardCode( const vec & xtheta, 
                                 const gpcov & CovV, 
                                 const gpcov & CovR, 
                                 const double & sigma, 
                                 const mat & yobs,
                                 const std::function<mat (vec, mat)> & fODE) {
  int n = (xtheta.size() - 3)/2;
  lp ret;
  ret.gradient.set_size(xtheta.size());
  
  const double *xthetaPtr = xtheta.memptr();
  const double *VmphiPtr = CovV.mphiBand.memptr();
  const double *VKinvPtr = CovV.KinvBand.memptr();
  const double *VCinvPtr = CovV.CinvBand.memptr();
  const double *RmphiPtr = CovR.mphiBand.memptr(); 
  const double *RKinvPtr = CovR.KinvBand.memptr(); 
  const double *RCinvPtr = CovR.CinvBand.memptr();
  const double *sigmaPtr = &sigma; 
  const double *yobsPtr = yobs.memptr();
  double *retPtr = &ret.value;
  double *retgradPtr = ret.gradient.memptr();
  
  xthetallikBandC( xthetaPtr, VmphiPtr, VKinvPtr, VCinvPtr,
                   RmphiPtr, RKinvPtr, RCinvPtr, &CovV.bandsize, &n,
                   sigmaPtr, yobsPtr, retPtr, retgradPtr, fODE);
  
  return ret;
}


//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
arma::vec speedbenchmarkXthetallik(const arma::mat & yobs, 
                                   const Rcpp::List & covVr, 
                                   const Rcpp::List & covRr, 
                                   const double & sigma, 
                                   const arma::vec & initial,
                                   const int & nrep = 10000){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta);
  
  std::vector<chrono::high_resolution_clock::time_point> timestamps;
  std::vector<lp> llikResults;
  
  // capture run time here
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallik_rescaled(initial, covV, covR, sigma, yobs, fnmodelODE));  
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallikBandApproxHardCode(initial, covV, covR, sigma, yobs, fnmodelODE));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallikHardCode(initial, covV, covR, sigma, yobs, fnmodelODE));  
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallik(initial, covV, covR, sigma, yobs, fnmodel));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallik_withmu(initial, covV, covR, sigma, yobs, fnmodel));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallik_withmu2(initial, covV, covR, sigma, yobs, fnmodel));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallikBandApprox(initial, covV, covR, sigma, yobs, fnmodel));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    llikResults.push_back(xthetallikWithmuBand(initial, covV, covR, sigma, yobs, fnmodel));
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  
  
  if(abs(llikResults[6].value - llikResults[1].value) > 1e-10
       || sum(abs(llikResults[6].gradient - llikResults[1].gradient)) > 1e-8){
    throw "xthetallikBandApprox and xthetallikBandApproxHardCode not agree";
  }
  
  arma::vec returnValues(timestamps.size()-1);
  for(int i = 0; i < timestamps.size()-1; i++){
    returnValues(i) = chrono::duration_cast<chrono::nanoseconds>(timestamps[i+1]-timestamps[i]).count();
  }
  return returnValues;
}
