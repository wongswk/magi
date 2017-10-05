#include <RcppArmadillo.h>

#include "hmc.h"
#include "classDefinition.h"

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