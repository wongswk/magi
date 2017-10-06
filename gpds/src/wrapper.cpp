#include <cmath>
#include <random>
#include <vector>
// [[Rcpp::plugins(cpp11)]]
#include <functional>

#include "classDefinition.h"
#include "hmc.h"
#include "tgtdistr.h"
#include "paralleltempering.h"

using namespace std;
using namespace arma;
using namespace Rcpp;


//' R wrapper for phisigllik
//' @export
// [[Rcpp::export]]
Rcpp::List phisigllikC(const arma::vec & phisig, 
                       const arma::mat & yobs, 
                       const arma::mat & dist, 
                       std::string kernel="matern"){
  lp ret = phisigllik(phisig, yobs, dist, kernel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' sample from GP marginal likelihood for phi and sigma
//' @export
// [[Rcpp::export]]
Rcpp::List phisigSample( const arma::mat & yobs, 
                         const arma::mat & dist, 
                         const arma::vec & initial, 
                         const arma::vec & step,
                         int nsteps = 1, 
                         bool traj = false, 
                         std::string kernel = "matern"){
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

gpcov cov_r2cpp(const Rcpp::List & cov_r){
  // TODO: gpcov shoule be a class and this should be a constructor
  // should only copy when needed
  // set member function and private member elements
  // initiate member elements only when get function called
  // store pointer in consructor
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
  cov_v.CinvBand = as<vec>(cov_r["CinvBand"]);
  cov_v.mphiBand = as<vec>(cov_r["mphiBand"]);
  cov_v.KinvBand = as<vec>(cov_r["KinvBand"]);
  cov_v.bandsize = as<int>(cov_r["bandsize"]);
  return cov_v;
}

//' sample from GP ODE for latent x and theta
//' @export
// [[Rcpp::export]]
Rcpp::List xthetaSample( const arma::mat & yobs, 
                         const Rcpp::List & covVr, 
                         const Rcpp::List & covRr, 
                         const double & sigma, 
                         const arma::vec & initial, 
                         const arma::vec & step,
                         const int nsteps = 1, 
                         const bool traj = false, 
                         const std::string loglikflag = "usual"){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  std::function<lp(vec)> tgt;
  if(loglikflag == "rescaled"){
    tgt = std::bind(xthetallik_rescaled, std::placeholders::_1, 
                    covV, covR, sigma, yobs, fnmodelODE);
  }else if(loglikflag == "usual"){
    tgt = std::bind(xthetallik, std::placeholders::_1, 
                    covV, covR, sigma, yobs, fnmodelODE);
  }else if(loglikflag == "band"){
    tgt = std::bind(xthetallikBandApprox, std::placeholders::_1, 
                    covV, covR, sigma, yobs);
  }else{
    throw "loglikflag must be 'usual', 'rescaled', or 'band'";
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

//' parallel tempered version of hmc for xtheta sample
//' @export
// [[Rcpp::export]]
arma::cube parallel_temper_hmc_xtheta( const arma::mat & yobs, 
                                       const Rcpp::List & covVr, 
                                       const Rcpp::List & covRr, 
                                       const double & sigma, 
                                       const arma::vec & temperature, 
                                       const double & alpha0, 
                                       const arma::vec & initial, 
                                       const arma::vec & step, 
                                       int nsteps = 1, 
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
    
  
  // cout << "test tgt value = " << tgt(initial).value 
  //      << " &tgt = " << &tgt << endl;
  // cout << "test & HMC func = " << &basic_hmcC << endl;
  
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
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallikC(const arma::mat & yobs, 
                       const Rcpp::List & covVr, 
                       const Rcpp::List & covRr, 
                       const double & sigma, 
                       const arma::vec & initial){
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  lp ret = xthetallik(initial, covV, covR, sigma, yobs, fnmodelODE);
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
  gpcov covV = cov_r2cpp(covVr);
  gpcov covR = cov_r2cpp(covRr);
  lp ret = xthetallikBandApprox(initial, covV, covR, sigma, yobs);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}
