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
#include "paralleltempering.h"
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

gpcov cov_r2cpp(List cov_r){
  gpcov cov_v;
  cov_v.C = as<mat>(cov_r["C"]);
  cov_v.Cinv = as<mat>(cov_r["Cinv"]);
  cov_v.mphi = as<mat>(cov_r["mphi"]);
  cov_v.Kphi = as<mat>(cov_r["Kphi"]);
  cov_v.Kinv = as<mat>(cov_r["Kinv"]);
  return cov_v;
}

//' sample from GP ODE for latent x and theta
// [[Rcpp::export]]
Rcpp::List xthetaSample( mat yobs, List covVr, List covRr, double sigma, const vec & initial, vec step,
                         int nsteps = 1, bool traj = false){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  std::function<lp(vec)> tgt = std::bind(xthetallik, std::placeholders::_1, 
                   covV, covR, sigma, yobs, fnmodelODE);
  vec lb = ones<vec>(initial.size()) * (-datum::inf);
  lb.subvec(lb.size() - 3, lb.size() - 1).fill(0.0);
  
  // cout << lb << endl;
  
  // lp tmp = tgt(initial);
  // cout << tmp.value << "\n" << tmp.gradient << endl;
  
  hmcstate post = basic_hmcC(tgt, initial, step, lb, {datum::inf}, nsteps, traj);
  
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


// [[Rcpp::export]]
arma::cube parallel_temper_hmc_xtheta( 
    mat yobs, List covVr, List covRr, double sigma, 
    const vec & initial, const vec & step, int nsteps = 1){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  std::function<lp(vec)> tgt = std::bind(xthetallik, std::placeholders::_1, 
                   covV, covR, sigma, yobs, fnmodelODE);
  
  vec lb = ones<vec>(initial.size()) * (-datum::inf);
  lb.subvec(lb.size() - 3, lb.size() - 1).fill(0.0);
  
  std::function<mcmcstate(std::function<lp(vec)>, mcmcstate)> hmc_simple =
    [&lb, &step, nsteps](std::function<lp(vec)> tgt_tempered, mcmcstate currstate){
      hmcstate post = basic_hmcC(tgt_tempered, currstate.state, step, lb, 
                                 {datum::inf}, nsteps, false);
      return mcmcstate(post);
    };
    
  vec temperature = {1, 1.5, 2, 3, 4.5, 6, 8};
  
  cube samples = parallel_termperingC(tgt, 
                                      hmc_simple, 
                                      temperature, 
                                      initial, 
                                      0.10, 
                                      1e5);
  return samples;
}

