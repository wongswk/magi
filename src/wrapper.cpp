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

using namespace std;
using namespace arma;
using namespace Rcpp;

//' R wrapper for basic_hmcC
// [[Rcpp::export]]
Rcpp::List hmc(const vec & initial, vec step,
               int nsteps = 1, bool traj = false){
  hmcstate post = basic_hmcC(lpnormal, initial, step, nsteps, traj);
  
  return List::create(Named("final")=post.final,
                      Named("final.p")=post.finalp,
                      Named("lpr")=post.lprvalue,
                      Named("step")=post.step,
                      Named("apr")=post.apr,
                      Named("acc")=post.acc,
                      Named("delta")=post.delta);
}

//' R wrapper for phisigllik
// [[Rcpp::export]]
Rcpp::List phisigllikTest(vec phisig, mat yobs, mat dist){
  lp ret = phisigllik(phisig, yobs, dist);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}
