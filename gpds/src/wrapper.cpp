// [[Rcpp::plugins(cpp11)]]
#include "wrapper.h"

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
  gpcov cov_v;
  
  // std::cout << "extract list as constant reference\n";
  const Rcpp::NumericMatrix & Cinv = as<const NumericMatrix>(cov_r["Cinv"]);
  const Rcpp::NumericMatrix & mphi = as<const NumericMatrix>(cov_r["mphi"]);
  const Rcpp::NumericMatrix & Kinv = as<const NumericMatrix>(cov_r["Kinv"]);
  const Rcpp::NumericMatrix & CinvBand = as<const NumericMatrix>(cov_r["CinvBand"]);
  const Rcpp::NumericMatrix & mphiBand = as<const NumericMatrix>(cov_r["mphiBand"]);
  const Rcpp::NumericMatrix & KinvBand = as<const NumericMatrix>(cov_r["KinvBand"]);
  const Rcpp::NumericVector & mu = as<const NumericVector>(cov_r["mu"]);
  const Rcpp::NumericVector & dotmu = as<const NumericVector>(cov_r["dotmu"]);
  
  // *(const_cast<double*>( &(mu[0]))) = -1; // this part is working -- R value changed
  
  // std::cout << "use R memory without copy\n";
  cov_v.Cinv = mat(const_cast<double*>( Cinv.begin()), Cinv.nrow(), Cinv.ncol(), false, false);
  cov_v.mphi = mat(const_cast<double*>( mphi.begin()), mphi.nrow(), mphi.ncol(), false, false);
  cov_v.Kinv = mat(const_cast<double*>( Kinv.begin()), Kinv.nrow(), Kinv.ncol(), false, false);
  cov_v.CinvBand = mat(const_cast<double*>( CinvBand.begin()), CinvBand.nrow(), CinvBand.ncol(), false, false);
  cov_v.mphiBand = mat(const_cast<double*>( mphiBand.begin()), mphiBand.nrow(), mphiBand.ncol(), false, false);
  cov_v.KinvBand = mat(const_cast<double*>( KinvBand.begin()), KinvBand.nrow(), KinvBand.ncol(), false, false);
  cov_v.mu = vec(const_cast<double*>( &(mu[0])), mu.size(), false, false);
  cov_v.dotmu = vec(const_cast<double*>( &(dotmu[0])), dotmu.size(), false, false);
  cov_v.bandsize = as<int>(cov_r["bandsize"]);
  
  // cov_v.mu(1) = 2; // this part is also working -- R value changed
  // cov_v.CinvBand(0) = 999;
  return cov_v;
}

//' sample from GP ODE for latent x and theta
//' @export
// [[Rcpp::export]]
Rcpp::List xthetaSample( const arma::mat & yobs, 
                         const Rcpp::List & covList, 
                         const double & sigma, 
                         const arma::vec & initial, 
                         const arma::vec & step,
                         const int nsteps = 1, 
                         const bool traj = false, 
                         const std::string loglikflag = "usual",
                         const double & overallTemperature = 1,
                         const double & priorTemperature = 1, 
                         const std::string modelName = "FN"){
  vector<gpcov> covAllDimensions(covList.size());
  for(unsigned int i = 0; i < covList.size(); i++){
    covAllDimensions[i] = cov_r2cpp(covList[i]);
  }
  std::function<lp(vec)> tgt;
  OdeSystem model;
  
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else{
    throw "modelName bust be one of 'FN', 'Hes1'";
  }
    
  if(loglikflag == "usual"){
    tgt = std::bind(xthetallik, std::placeholders::_1, 
                    covAllDimensions, sigma, yobs, model, false, priorTemperature);
  }else if(loglikflag == "withmean"){
    tgt = std::bind(xthetallikWithmuBand, std::placeholders::_1, 
                    covAllDimensions, sigma, yobs, model, false, priorTemperature);
  }else if(loglikflag == "band"){
    tgt = std::bind(xthetallik, std::placeholders::_1, 
                    covAllDimensions, sigma, yobs, model, true, priorTemperature);
  }else if(loglikflag == "withmeanBand"){
    tgt = std::bind(xthetallikWithmuBand, std::placeholders::_1, 
                    covAllDimensions, sigma, yobs, model, true, priorTemperature);
  }else{
    throw "loglikflag must be 'usual', 'withmean', 'band', or 'withmeanBand'";
  }
  vec lb = ones<vec>(initial.size()) * (-datum::inf);
  lb.subvec(yobs.size(), lb.size() - 1) = model.thetaLowerBound;
  
  // cout << lb << endl;
  
  // lp tmp = tgt(initial);
  // cout << tmp.value << "\n" << tmp.gradient << endl;
  std::function<lp(vec)> tgtTempered = [&tgt, &overallTemperature](vec xInput) -> lp {
    lp ret = tgt(xInput);
    ret.value /= overallTemperature;
    ret.gradient /= overallTemperature;
    return ret;
  };
  
  hmcstate post = basic_hmcC(tgtTempered, initial, step, lb, {datum::inf}, nsteps, traj);
  
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
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  
  
  std::function<lp(vec)> tgt = std::bind(xthetallik, std::placeholders::_1, 
                   covAllDimensions, sigma, yobs, fnmodel, true, 1.0);
  
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
  // cout << "test hmc_simple = " << initpost_mcmcstate.lpv << endl
  //      << initpost_mcmcstate.acc << endl;
       // << initpost_mcmcstate.state << endl;
  
  // cout << "prepare to call parallel_termperingC" << endl;
  
  cube samples = parallel_termperingC(tgt,
                                      hmc_simple,
                                      temperature,
                                      initial,
                                      alpha0,
                                      niter,
                                      true);
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
                       const arma::vec & initial,
                       const bool useBand = false,
                       const double & priorTemperature = 1.0){
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  lp ret = xthetallik(initial, covAllDimensions, sigma, yobs, fnmodel, useBand, priorTemperature);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}


//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallikWithmuBandC(const arma::mat & yobs, 
                                 const Rcpp::List & covVr, 
                                 const Rcpp::List & covRr, 
                                 const double & sigma, 
                                 const arma::vec & initial,
                                 const bool useBand = true,
                                 const double & priorTemperature = 1.0){
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  lp ret = xthetallikWithmuBand(initial, covAllDimensions, sigma, yobs, fnmodel, useBand, priorTemperature);
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}
