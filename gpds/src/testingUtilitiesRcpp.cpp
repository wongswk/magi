#include <RcppArmadillo.h>
#include <chrono>

#include "hmc.h"
#include "classDefinition.h"
#include "wrapper.h"
#include "dynamicalSystemModels.h"
#include "band.h"

using namespace arma;
using namespace Rcpp;
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

gpcov cov_r2cpp_legacy(const Rcpp::List & cov_r){
  gpcov cov_v;
  cov_v.C = as<mat>(cov_r["C"]);
  cov_v.Cinv = as<mat>(cov_r["Cinv"]);
  cov_v.mphi = as<mat>(cov_r["mphi"]);
  cov_v.Kphi = as<mat>(cov_r["Kphi"]);
  cov_v.Kinv = as<mat>(cov_r["Kinv"]);
  cov_v.Ceigen1over = as<vec>(cov_r["Ceigen1over"]);
  cov_v.Keigen1over = as<vec>(cov_r["Keigen1over"]);
  cov_v.CeigenVec = as<mat>(cov_r["CeigenVec"]);
  cov_v.KeigenVec = as<mat>(cov_r["KeigenVec"]);
  cov_v.CinvBand = as<mat>(cov_r["CinvBand"]);
  cov_v.mphiBand = as<mat>(cov_r["mphiBand"]);
  cov_v.KinvBand = as<mat>(cov_r["KinvBand"]);
  cov_v.mu = as<vec>(cov_r["mu"]);
  cov_v.dotmu = as<vec>(cov_r["dotmu"]);
  cov_v.bandsize = as<int>(cov_r["bandsize"]);
  return cov_v;
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

lp xthetallikTwoDimension( const vec & xtheta, 
                           const gpcov & CovV, 
                           const gpcov & CovR, 
                           const double & sigma, 
                           const mat & yobs, 
                           const OdeSystem & fOdeModel) {
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
  
  const mat & fderiv = fOdeModel.fOde(theta, join_horiz(Vsm, Rsm));
  const cube & fderivDx = fOdeModel.fOdeDx(theta, join_horiz(Vsm, Rsm));
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, join_horiz(Vsm, Rsm));
  
  mat res(2,3);
  
  // V 
  vec frV = (fderiv.col(0) - CovV.mphi * Vsm); // n^2 operation
  vec fitLevelErrorV = Vsm - yobs.col(0);
  fitLevelErrorV(find_nonfinite(fitLevelErrorV)).fill(0.0);
  res(0,0) = -0.5 * sum(square( fitLevelErrorV )) / pow(sigma,2);
  res(0,1) = -0.5 * as_scalar( frV.t() * CovV.Kinv * frV);
  res(0,2) = -0.5 * as_scalar( Vsm.t() * CovV.Cinv * Vsm);
  
  // R
  vec frR = (fderiv.col(1) - CovR.mphi * Rsm); // n^2 operation
  vec fitLevelErrorR = Rsm - yobs.col(1);
  fitLevelErrorR(find_nonfinite(fitLevelErrorR)).fill(0.0);
  
  res(1,0) = -0.5 * sum(square( fitLevelErrorR )) / pow(sigma,2);
  res(1,1) = -0.5 * as_scalar( frR.t() * CovR.Kinv * frR);
  res(1,2) = -0.5 * as_scalar( Rsm.t() * CovR.Cinv * Rsm);
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = -CovV.mphi;
  Vtemp.diag() += fderivDx.slice(0).col(0);
  
  vec KinvFrV = (CovV.Kinv * frV);
  
  vec VC2 =  2.0 * join_vert(join_vert( Vtemp.t() * KinvFrV, // n^2 operation
                                        fderivDx.slice(0).col(1) % KinvFrV ),
                                        fderivDtheta.slice(0).t() * KinvFrV );
  
  
  // R contrib
  mat Rtemp = -CovR.mphi;
  Rtemp.diag() += fderivDx.slice(1).col(1);
  
  vec KinvFrR = (CovR.Kinv * frR);
  vec RC2 = 2.0 * join_vert(join_vert( fderivDx.slice(1).col(0) % KinvFrR,
                                       Rtemp.t() * KinvFrR), // n^2 operation
                                       fderivDtheta.slice(1).t() * KinvFrR );
  
  vec C3 = join_vert(join_vert( 2.0 * CovV.Cinv * Vsm,  
                                2.0 * CovR.Cinv * Rsm ), 
                                zeros<vec>(theta.size()));  
  vec C1 = join_vert(join_vert( 2.0 * fitLevelErrorV / pow(sigma,2) ,  
                                2.0 * fitLevelErrorR / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
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
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp_legacy(covVr);
  covAllDimensions[1] = cov_r2cpp_legacy(covRr);
  
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  
  std::vector<chrono::high_resolution_clock::time_point> timestamps;
  
  int nrepShort = nrep/100;
  
  // capture run time here
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrepShort; i++){
    lp tmp1 = xthetallik_rescaled(initial, covAllDimensions[0], covAllDimensions[1], sigma, yobs, fnmodelODE);  
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp2 = xthetallikBandApproxHardCode(initial, covAllDimensions[0], covAllDimensions[1], sigma, yobs, fnmodelODE);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp3 = xthetallikHardCode(initial, covAllDimensions[0], covAllDimensions[1], sigma, yobs, fnmodelODE);  
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp4 = xthetallik(initial, covAllDimensions, sigma, yobs, fnmodel);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp5 = xthetallik_withmu(initial, covAllDimensions, sigma, yobs, fnmodel);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp6 = xthetallikWithmuBand(initial, covAllDimensions, sigma, yobs, fnmodel, false);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp7 = xthetallikBandApprox(initial, covAllDimensions, sigma, yobs, fnmodel);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp8 = xthetallikWithmuBand(initial, covAllDimensions, sigma, yobs, fnmodel, true);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  for(int i=0; i < nrep; i++){
    lp tmp9 = xthetallikTwoDimension(initial, covAllDimensions[0], covAllDimensions[1], sigma, yobs, fnmodel);
  }
  timestamps.push_back(chrono::high_resolution_clock::now());
  
  arma::vec returnValues(timestamps.size()-1);
  for(unsigned int i = 0; i < timestamps.size()-1; i++){
    returnValues(i) = chrono::duration_cast<chrono::nanoseconds>(timestamps[i+1]-timestamps[i]).count();
  }
  returnValues /= nrep;
  returnValues(0) *= nrep/nrepShort;
  return returnValues;
}

//' test for pass by reference for gpcov
//' @export
// [[Rcpp::export]]
int changeGPcovFromC(Rcpp::List & covVr){
  gpcov covV = cov_r2cpp(covVr);
  covV.Cinv.fill(1);
  covV.mphi.fill(2);
  covV.Kinv.fill(3);
  covV.CinvBand.fill(4);
  covV.mphiBand.fill(5);
  covV.KinvBand.fill(6);
  covV.mu.fill(77);
  covV.dotmu.fill(666);
  return 0;
}

// [[Rcpp::export]]
void cov_r2cpp_t1(const Rcpp::List & cov_r){
  const double* CinvPtr = as<const NumericMatrix>(cov_r["Cinv"]).begin();
  *(const_cast<double*> (CinvPtr) )= 0;
  // std::cout << CinvPtr << endl;
  // std::cout << as<mat>(cov_r["Cinv"]).memptr() << endl;
}

// [[Rcpp::export]]
void cov_r2cpp_t2(Rcpp::NumericMatrix & cov_r){
  // std::cout << cov_r.begin() << endl;
  // std::cout << as<const NumericMatrix>(cov_r).begin() << endl;
  // std::cout << as<mat>(cov_r).memptr();
  cov_r[0] = 0;
}

// [[Rcpp::export]]
void cov_r2cpp_t3(arma::mat & cov_r){
  // std::cout << cov_r.memptr() << endl;
  cov_r(0) = 0;
}

