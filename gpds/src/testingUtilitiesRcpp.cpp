#include <RcppArmadillo.h>
#include <chrono>

#include "hmc.h"
#include "classDefinition.h"
#include "wrapper.h"

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
  // capture run time here
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  for(int i=0; i < nrep; i++){
    lp ret1 = xthetallik_rescaled(initial, covV, covR, sigma, yobs, fnmodelODE);  
  }
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  for(int i=0; i < nrep; i++){
    lp ret2 = xthetallikBandApprox(initial, covV, covR, sigma, yobs);
  }
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  for(int i=0; i < nrep; i++){
    lp ret3 = xthetallik(initial, covV, covR, sigma, yobs, fnmodelODE);  
  }
  chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
  double duration1 = chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();
  double duration2 = chrono::duration_cast<chrono::nanoseconds>(t2-t1).count();
  double duration3 = chrono::duration_cast<chrono::nanoseconds>(t3-t2).count();
  return arma::vec({duration1, duration2, duration3});
}
