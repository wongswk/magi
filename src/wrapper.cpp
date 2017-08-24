// #include <Rcpp.h>
// using namespace Rcpp;

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include "hmc.h"
#include "tgtdistr.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
#include <functional>

using namespace std;
using namespace arma;
using namespace Rcpp;

//' R wrapper for basic_hmcC
// [[Rcpp::export]]
Rcpp::List hmc(const vec & initial, vec step, vec lb, vec ub,
               int nsteps = 1, bool traj = false){
  // cout << lb << "\t" << ub << endl;
  hmcstate post = basic_hmcC(lpnormal, initial, step, 
                             lb, 
                             ub, 
                             nsteps, traj);
  Rcpp::List ret = List::create(Named("final")=post.final,
                                Named("final.p")=post.finalp,
                                Named("lpr")=post.lprvalue,
                                Named("step")=post.step,
                                Named("apr")=post.apr,
                                Named("acc")=post.acc,
                                Named("delta")=post.delta);
  if(traj){
    ret.push_back(post.trajp, "traj.p");
    ret.push_back(post.trajq, "traj.q");
    ret.push_back(post.trajH, "traj.H");
  }
  return ret;
}

//' R wrapper for phisigllik
// [[Rcpp::export]]
Rcpp::List phisigllikTest(vec phisig, mat yobs, mat dist){
  lp ret = phisigllik(phisig, yobs, dist);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' sample from GP marginal likelihood for phi and sigma
// [[Rcpp::export]]
Rcpp::List phisigSample( mat yobs, mat dist, const vec & initial, vec step,
                         int nsteps = 1, bool traj = false){
  std::function<lp(vec)> tgt = std::bind(phisigllik, std::placeholders::_1, yobs, dist);
  hmcstate post = basic_hmcC(tgt, initial, step, 
                             std::vector<double>({0.0}), 
                             std::vector<double>({datum::inf}), 
                             nsteps, traj);
  lp ret = tgt(initial);
  return List::create(Named("final")=post.final,
                      Named("final.p")=post.finalp,
                      Named("lpr")=post.lprvalue,
                      Named("step")=post.step,
                      Named("apr")=post.apr,
                      Named("acc")=post.acc,
                      Named("delta")=post.delta);
}

