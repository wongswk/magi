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
Rcpp::List phisigllikTest(vec phisig, mat yobs, mat dist, string kernel="matern"){
  lp ret = phisigllik(phisig, yobs, dist, kernel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' sample from GP marginal likelihood for phi and sigma
// [[Rcpp::export]]
Rcpp::List phisigSample( mat yobs, mat dist, const vec & initial, vec step,
                         int nsteps = 1, bool traj = false, string kernel = "matern"){
  std::function<lp(vec)> tgt = std::bind(phisigllik, std::placeholders::_1, yobs, dist, kernel);
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
  cov_v.Ceigen1over = as<vec>(cov_r["Ceigen1over"]);
  cov_v.Keigen1over = as<vec>(cov_r["Keigen1over"]);
  cov_v.CeigenVec = as<mat>(cov_r["CeigenVec"]);
  cov_v.KeigenVec = as<mat>(cov_r["KeigenVec"]);
  return cov_v;
}

//' sample from GP ODE for latent x and theta
// [[Rcpp::export]]
Rcpp::List xthetaSample( const mat & yobs, const List & covVr, const List & covRr, double sigma, const vec & initial, vec step,
                         int nsteps = 1, bool traj = false, bool rescaleloglik = false){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  std::function<lp(vec)> tgt;
  if(rescaleloglik){
    tgt = std::bind(xthetallik_rescaled, std::placeholders::_1, 
                    covV, covR, sigma, yobs, fnmodelODE);
  }else{
    tgt = std::bind(xthetallik, std::placeholders::_1, 
                    covV, covR, sigma, yobs, fnmodelODE);
  }
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
    mat yobs, List covVr, List covRr, double sigma, const vec & temperature, 
    const double & alpha0, const vec & initial, const vec & step, int nsteps = 1, 
    int niter=1e4){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  std::function<lp(vec)> tgt = std::bind(xthetallik, std::placeholders::_1, 
                   covV, covR, sigma, yobs, fnmodelODE);
  
  vec lb = ones<vec>(initial.size()) * (-datum::inf);
  lb.subvec(lb.size() - 3, lb.size() - 1).fill(0.0);
  
  std::function<mcmcstate(std::function<lp(vec)>, mcmcstate)> hmc_simple =
    [step, lb, nsteps](std::function<lp(vec)> tgt_tempered, mcmcstate currstate) -> mcmcstate{
      currstate.lpv = tgt_tempered(currstate.state).value;
      vec rstep = arma::randu<vec>(step.size()) % step  + step;
      hmcstate post = basic_hmcC(tgt_tempered, currstate.state, rstep, lb,
                                 {arma::datum::inf}, nsteps, true);
      // cout << post.trajH;
      return mcmcstate(post);
    };
    
  
  cout << "test tgt value = " << tgt(initial).value 
       << " &tgt = " << &tgt << endl;
  cout << "test & HMC func = " << &basic_hmcC << endl;
  
  mcmcstate init_mcmcstate;
  init_mcmcstate.state = initial;
  mcmcstate initpost_mcmcstate = hmc_simple(tgt, init_mcmcstate);
  cout << "test hmc_simple = " << initpost_mcmcstate.lpv << endl
       << initpost_mcmcstate.acc << endl;
       // << initpost_mcmcstate.state << endl;
  
  cout << "prepare to call parallel_termperingC" << endl;
  
  cube samples = parallel_termperingC(tgt,
                                      hmc_simple,
                                      temperature,
                                      initial,
                                      alpha0,
                                      niter);
  return samples;
  // return arma::zeros<cube>(1,1,1);
}

//' R wrapper for xthetallik
// [[Rcpp::export]]
Rcpp::List xthetallikTest(mat yobs, List covVr, List covRr, double sigma, const vec & initial){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  lp ret = xthetallik(initial, covV, covR, sigma, yobs, fnmodelODE);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}
