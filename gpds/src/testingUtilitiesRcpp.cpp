#include <RcppArmadillo.h>
#include <chrono>

#include "hmc.h"
#include "classDefinition.h"
#include "wrapper.h"
#include "dynamicalSystemModels.h"
#include "band.h"

using namespace arma;
using namespace Rcpp;

lp xthetallik_withmu( const arma::vec & xtheta, 
                      const std::vector<gpcov> & CovAllDimensions,
                      const double & sigma, 
                      const arma::mat & yobs, 
                      const OdeSystem & fOdeModel);
lp xthetallik_rescaled( const arma::vec & xtheta,
                        const gpcov & CovV,
                        const gpcov & CovR,
                        const double & sigma,
                        const arma::mat & yobs,
                        const std::function<arma::mat (arma::vec, arma::mat)> & fODE);
lp xthetallikBandApprox( const arma::vec & xtheta, 
                         const std::vector<gpcov> & CovAllDimensions,
                         const double & sigma, 
                         const arma::mat & yobs,
                         const OdeSystem & fOdeModel);
lp xthetallikTwoDimension( const vec & xtheta, 
                           const gpcov & CovV, 
                           const gpcov & CovR, 
                           const double & sigma, 
                           const mat & yobs, 
                           const OdeSystem & fOdeModel);
lp xthetallikBandApproxHardCode( const vec & xtheta, 
                                 const gpcov & CovV, 
                                 const gpcov & CovR, 
                                 const double & sigma, 
                                 const mat & yobs,
                                 const std::function<mat (vec, mat)> & fODE);
lp xthetallikHardCode( const vec & xtheta, 
                       const gpcov & CovV, 
                       const gpcov & CovR, 
                       const double & sigma, 
                       const mat & yobs, 
                       const std::function<mat (vec, mat)> & fODE);

lp phisigllikHard2D( const arma::vec &, const arma::mat &, const arma::mat &, string kernel = "matern");

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




//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallik_rescaledC(const arma::mat & yobs, 
                                const Rcpp::List & covVr, 
                                const Rcpp::List & covRr, 
                                const double & sigma, 
                                const arma::vec & initial){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  lp ret = xthetallik_rescaled(initial, covV, covR, sigma, yobs, fnmodelODE);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' R wrapper for xthetallikBandApprox
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallikBandApproxC( arma::mat & yobs, 
                                  const Rcpp::List & covVr, 
                                  const Rcpp::List & covRr, 
                                  double & sigma, 
                                  arma::vec & initial){
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  lp ret = xthetallikBandApprox(initial, covAllDimensions, sigma, yobs, fnmodel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallik_withmuC(const arma::mat & yobs, 
                              const Rcpp::List & covVr, 
                              const Rcpp::List & covRr, 
                              const double & sigma, 
                              const arma::vec & initial){
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  lp ret = xthetallik_withmu(initial, covAllDimensions, sigma, yobs, fnmodel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
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
    lp tmp4 = xthetallik(initial, covAllDimensions, sigma, yobs, fnmodel, false);
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
  for(int i=0; i < nrep; i++){
    lp tmp10 = xthetallik(initial, covAllDimensions, sigma, yobs, fnmodel, true);
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

//' R wrapper for phisigllik
//' @export
// [[Rcpp::export]]
Rcpp::List phisigllikHard2DC(const arma::vec & phisig, 
                             const arma::mat & yobs, 
                             const arma::mat & dist, 
                             std::string kernel="matern"){
  lp ret = phisigllikHard2D(phisig, yobs, dist, kernel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}
