// [[Rcpp::plugins(cpp11)]]
#include "RcppWrapper.h"
#include "RcppArmadillo.h"
#include "Sampler.h"
#include "gpsmoothing.h"

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

//' R wrapper for phisigllik
//' @export
// [[Rcpp::export]]
Rcpp::List phisigloocvllikC(const arma::vec & phisig, 
                            const arma::mat & yobs, 
                            const arma::mat & dist, 
                            std::string kernel="matern"){
  lp ret = phisigloocvllik(phisig, yobs, dist, kernel);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}

//' R wrapper for phisigloocvmse
//' @export
// [[Rcpp::export]]
Rcpp::List phisigloocvmseC(const arma::vec & phisig, 
                           const arma::mat & yobs, 
                           const arma::mat & dist, 
                           std::string kernel="matern"){
  lp ret = phisigloocvmse(phisig, yobs, dist, kernel);
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
  const Rcpp::NumericVector & tvecCovInput = as<const NumericVector>(cov_r["tvecCovInput"]);
  
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
  cov_v.tvecCovInput = vec(const_cast<double*>( &(tvecCovInput[0])), tvecCovInput.size(), false, false);
  
  // cov_v.mu(1) = 2; // this part is also working -- R value changed
  // cov_v.CinvBand(0) = 999;
  return cov_v;
}


Rcpp::List cov_cpp2r(const gpcov & cov_v){
  // TODO: gpcov shoule be a class and this should be a constructor
  Rcpp::List cov_r = List::create(
    Named("C")=cov_v.C,
    Named("dCdphiCube")=cov_v.dCdphiCube,
    Named("Cprime")=cov_v.Cprime,
    Named("Cdoubleprime")=cov_v.Cdoubleprime,
    Named("dCprimedphiCube")=cov_v.dCprimedphiCube,
    Named("dCdoubleprimedphiCube")=cov_v.dCdoubleprimedphiCube,
    Named("Cinv")=cov_v.Cinv,
    Named("mphi")=cov_v.mphi,
    Named("Kinv")=cov_v.Kinv,
    Named("Sigma")=cov_v.Sigma,
    Named("dSigmadphiCube")=cov_v.dSigmadphiCube,
    Named("tvecCovInput")=cov_v.tvecCovInput
  );
  return cov_r;
}

//' general matern cov calculation Rcpp wrapper
//' @export
// [[Rcpp::export]]
Rcpp::List generalMaternCovRcpp( const arma::vec & phi, const arma::mat & distSigned, int complexity = 3){
  return cov_cpp2r(generalMaternCov(phi, distSigned, complexity));
}

//' sample from GP ODE for latent x and theta
//' @export
// [[Rcpp::export]]
Rcpp::List xthetaSample( const arma::mat & yobs, 
                         const Rcpp::List & covList, 
                         const arma::vec & sigmaInput,
                         const arma::vec & initial, 
                         const arma::vec & step,
                         const int nsteps = 1, 
                         const bool traj = false, 
                         const std::string loglikflag = "usual",
                         const double & overallTemperature = 1,
                         const Rcpp::NumericVector & priorTemperatureInput = 1.0, 
                         const std::string modelName = "FN"){
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  vec sigma( yobs.n_cols);
  if(sigmaInput.size() == 1){
    sigma.fill( as_scalar( sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput size not right");
  }
  
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
  }else if(modelName == "Hes1-log"){
    model = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
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
    throw std::runtime_error("loglikflag must be 'usual', 'withmean', 'band', or 'withmeanBand'");
  }
  vec lb = ones<vec>(initial.size()) * (-datum::inf);
  lb.subvec(yobs.size(), lb.size() - 1) = model.thetaLowerBound;
  
  // cout << lb << endl;
  
  // lp tmp = tgt(initial);
  // cout << tmp.value << "\n" << tmp.gradient << endl;
  std::function<lp(vec)> tgtTempered = [&tgt, &overallTemperature](const vec & xInput) -> lp {
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
//'
//' not used in the final method, in final method, only one temperature with heating is needed
//'
//' @noRd
// DISABLED [[Rcpp::export]]
arma::cube parallel_temper_hmc_xtheta( const arma::mat & yobs, 
                                       const Rcpp::List & covVr, 
                                       const Rcpp::List & covRr, 
                                       const arma::vec & sigmaInput,
                                       const arma::vec & temperature, 
                                       const double & alpha0, 
                                       const arma::vec & initial, 
                                       const arma::vec & step, 
                                       int nsteps = 1, 
                                       int niter=1e4){
  vec sigma( yobs.n_cols);
  if(sigmaInput.size() == 1){
    sigma.fill( as_scalar( sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput size not right");
  }
  
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  
  
  std::function<lp(vec)> tgt = std::bind(xthetallik, std::placeholders::_1, 
                   covAllDimensions, sigma, yobs, fnmodel, true, arma::ones(1));
  
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
Rcpp::List xthetallikRcpp(const arma::mat & yobs, 
                          const Rcpp::List & covAllDimInput,
                          const arma::vec & sigmaInput,
                          const arma::vec & initial,
                          const std::string modelName = "FN",
                          const bool useBand = false,
                          const Rcpp::NumericVector & priorTemperatureInput = 1.0){
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  vec sigma( yobs.n_cols);
  if(sigmaInput.size() == 1){
    sigma.fill( as_scalar( sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput size not right");
  }
  
  vector<gpcov> covAllDimensions(yobs.n_cols);
  for(unsigned j = 0; j < yobs.n_cols; j++){
    covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
  }
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "Hes1-log"){
    model = OdeSystem(hes1logmodelODE, hes1logmodelDx, hes1logmodelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'Hes1-log', 'HIV'");
  }
  
  lp ret = xthetallik(initial, covAllDimensions, sigma, yobs, model, useBand, priorTemperature);
  return List::create(Named("value")=ret.value,
                      Named("grad")=ret.gradient);
}


//' R wrapper for xthetallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetallikWithmuBandC(const arma::mat & yobs, 
                                 const Rcpp::List & covVr, 
                                 const Rcpp::List & covRr, 
                                 const arma::vec & sigmaInput,
                                 const arma::vec & initial,
                                 const bool useBand = true,
                                 const Rcpp::NumericVector & priorTemperatureInput = 1.0){
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  vec sigma( yobs.n_cols);
  if(sigmaInput.size() == 1){
    sigma.fill( as_scalar( sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput size not right");
  }
  vector<gpcov> covAllDimensions(2);
  covAllDimensions[0] = cov_r2cpp(covVr);
  covAllDimensions[1] = cov_r2cpp(covRr);
  OdeSystem fnmodel(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  lp ret = xthetallikWithmuBand(initial, covAllDimensions, sigma, yobs, fnmodel, useBand, priorTemperature);
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}


lp lp_r2cpp(const Rcpp::List & lp_r){
    lp lp_c;
    lp_c.value = lp_r["value"];

    const Rcpp::NumericVector & gradient = as<const NumericVector>(lp_r["gradient"]);
    lp_c.gradient = vec(const_cast<double*>( &(gradient[0])), gradient.size(), false, false);

    return lp_c;
}


//' R wrapper for basic_hmcC
//' @export
// [[Rcpp::export]]
Rcpp::List basic_hmcRcpp(const Rcpp::Function rlpr,
                         const arma::vec & initial,
                         const arma::vec & step,
                         arma::vec lb,
                         arma::vec ub,
                         const int nsteps = 1,
                         const bool traj = false){
    const std::function<lp (vec)> & clpr = [& rlpr](const vec & x) -> lp {
        return lp_r2cpp(rlpr(x));
    };
    const hmcstate & ret = basic_hmcC(clpr, initial, step, lb, ub, nsteps, traj);
    return Rcpp::List::create(Rcpp::Named("final")=ret.final,
                              Rcpp::Named("finalp")=ret.finalp,
                              Rcpp::Named("step")=ret.step,
                              Rcpp::Named("trajH")=ret.trajH,
                              Rcpp::Named("lprvalue")=ret.lprvalue,
                              Rcpp::Named("apr")=ret.apr,
                              Rcpp::Named("delta")=ret.delta,
                              Rcpp::Named("acc")=ret.acc,
                              Rcpp::Named("trajq")=ret.trajq,
                              Rcpp::Named("trajp")=ret.trajp);
}

//' R wrapper for chainSamplerRcpp
//' @export
// [[Rcpp::export]]
Rcpp::List chainSamplerRcpp(const arma::mat & yobs,
                            const Rcpp::List & covAllDimInput,
                            const int nstepsInput,
                            const std::string loglikflagInput,
                            const arma::vec & priorTemperatureInput,
                            const unsigned int sigmaSizeInput,
                            const OdeSystem & modelInput,
                            const unsigned int niterInput,
                            const double burninRatioInput,
                            const arma::vec &xthetasigmaInit,
                            const arma::vec &stepLowInit,
                            bool verbose){
  vector<gpcov> covAllDimensions(yobs.n_cols);
  for(unsigned j = 0; j < yobs.n_cols; j++){
    covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
  }
  Sampler sampler(yobs,
                  covAllDimensions,
                  nstepsInput,
                  loglikflagInput,
                  priorTemperatureInput,
                  sigmaSizeInput,
                  modelInput,
                  niterInput,
                  burninRatioInput);
  sampler.sampleChian(xthetasigmaInit, stepLowInit, verbose);
  return List::create(
            Named("lliklist")=sampler.lliklist,
            Named("xth")=sampler.xth
    );
}

//' R wrapper for chainSamplerRcpp
//' @export
// [[Rcpp::export]]
arma::vec optimizeThetaInitRcpp(const arma::mat & yobs,
                                const OdeSystem & modelInput,
                                const Rcpp::List & covAllDimInput,
                                const arma::vec & sigmaAllDimensionsInput,
                                const arma::vec & priorTemperatureInput,
                                const arma::mat & xInitInput,
                                const bool useBandInput){
    vector<gpcov> covAllDimensions(yobs.n_cols);
    for(unsigned j = 0; j < yobs.n_cols; j++){
        covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
    }
    return optimizeThetaInit(yobs,
                             modelInput,
                             covAllDimensions,
                             sigmaAllDimensionsInput,
                             priorTemperatureInput,
                             xInitInput,
                             useBandInput);
}
