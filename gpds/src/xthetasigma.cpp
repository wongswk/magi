#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
#include "wrapper.h"
#include "hmc.h"
#include "band.h"
#include "dynamicalSystemModels.h"

using namespace arma;


//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetasigmallik( const mat & xlatent, 
                    const vec & theta, 
                    const vec & sigmaInput, 
                    const mat & yobs, 
                    const std::vector<gpcov> & CovAllDimensions,
                    const OdeSystem & fOdeModel,
                    const arma::vec & priorTemperatureInput = ones(1),
                    const bool useBand = false,
                    const bool useMean = false) {
  
  if(useMean){
    mat xlatentShifted = xlatent;
    mat yobsShifted = yobs;
    mat muAllDimension(yobs.n_rows, yobs.n_cols);
    mat dotmuAllDimension(yobs.n_rows, yobs.n_cols);

    for(unsigned int i = 0; i < yobs.n_cols; i++){
      xlatentShifted.col(i) -= CovAllDimensions[i].mu;
      yobsShifted.col(i) -= CovAllDimensions[i].mu;
      muAllDimension.col(i) = CovAllDimensions[i].mu;
      dotmuAllDimension.col(i) = CovAllDimensions[i].dotmu;
    }
    
    OdeSystem fOdeModelShifted = fOdeModel;

    fOdeModelShifted.fOde = [&muAllDimension, &dotmuAllDimension, &fOdeModel]
    (const vec & theta, const mat & x) -> mat{
      return fOdeModel.fOde(theta, x+muAllDimension) - dotmuAllDimension;
    };
    
    fOdeModelShifted.fOdeDx = [&muAllDimension, &dotmuAllDimension, &fOdeModel]
    (const vec & theta, const mat & x) -> cube{ 
      return fOdeModel.fOdeDx(theta, x+muAllDimension);
    };
    
    fOdeModelShifted.fOdeDtheta = [&muAllDimension, &dotmuAllDimension, &fOdeModel]
    (const vec & theta, const mat & x) -> cube{ 
      return fOdeModel.fOdeDtheta(theta, x+muAllDimension);
    };

    return xthetasigmallik(xlatentShifted, theta, sigmaInput, yobsShifted, CovAllDimensions, 
                           fOdeModelShifted, priorTemperatureInput, useBand, false); 
  }
  
  int n = yobs.n_rows;
  int pdimension = yobs.n_cols;
  
  lp ret;
  if (fOdeModel.checkBound(xlatent, theta, &ret)) {
    return ret;
  }
  
  arma::vec priorTemperature(2);
  if(priorTemperatureInput.n_rows == 1){
    priorTemperature.fill(as_scalar(priorTemperatureInput));
  }else if(priorTemperatureInput.n_rows == 2){
    priorTemperature = priorTemperatureInput;
  }else{
    throw std::invalid_argument("priorTemperatureInput must be scaler or 2-vector");
  }
  
  vec sigma(yobs.n_cols);
  bool sigmaIsScaler = (sigmaInput.size() == 1);
  if (sigmaIsScaler){
    sigma.fill(as_scalar(sigmaInput));
  }else if(sigmaInput.size() == yobs.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput dimension not right");
  }
  vec sigmaSq = square(sigma);
  
  const mat & fderiv = fOdeModel.fOde(theta, xlatent);
  const cube & fderivDx = fOdeModel.fOdeDx(theta, xlatent);
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, xlatent);
  
  mat res(pdimension, 3);
  
  // V 
  mat fitLevelError = xlatent - yobs;
  mat fitDerivError(n, pdimension);
  vec nobs(pdimension);
  for( int vEachDim = 0; vEachDim < pdimension; vEachDim++){
    if(useBand){
      bmatvecmult(CovAllDimensions[vEachDim].mphiBand.memptr(), 
                  xlatent.colptr(vEachDim), 
                  &(CovAllDimensions[vEachDim].bandsize), 
                  &n, 
                  fitDerivError.colptr(vEachDim));
      fitDerivError.col(vEachDim) = fderiv.col(vEachDim) - fitDerivError.col(vEachDim);
    }else{
      fitDerivError.col(vEachDim) = fderiv.col(vEachDim);
      fitDerivError.col(vEachDim) -= CovAllDimensions[vEachDim].mphi * xlatent.col(vEachDim); // n^2 operation  
    }
    nobs(vEachDim) = uvec(find_finite(yobs.col(vEachDim))).size();
  }
  
  fitLevelError(find_nonfinite(fitLevelError)).fill(0.0);
  res.col(0) = -0.5 * sum(square( fitLevelError )).t() / sigmaSq - log(sigma) % nobs;
  vec sigmaGradient = sum(square( fitLevelError )).t() / (sigmaSq % sigma) - nobs / sigma;
  
  mat KinvfitDerivError(n, pdimension);
  mat CinvX(n, pdimension);
  
  for( int vEachDim = 0; vEachDim < pdimension; vEachDim++){
    if(useBand){
      bmatvecmult(CovAllDimensions[vEachDim].KinvBand.memptr(), 
                  fitDerivError.colptr(vEachDim), 
                  &(CovAllDimensions[vEachDim].bandsize), 
                  &n, 
                  KinvfitDerivError.colptr(vEachDim));
      bmatvecmult(CovAllDimensions[vEachDim].CinvBand.memptr(),
                  xlatent.colptr(vEachDim),  
                  &(CovAllDimensions[vEachDim].bandsize), 
                  &n, 
                  CinvX.colptr(vEachDim));
    }else{
      KinvfitDerivError.col(vEachDim) = CovAllDimensions[vEachDim].Kinv * fitDerivError.col(vEachDim);
      CinvX.col(vEachDim) = CovAllDimensions[vEachDim].Cinv * xlatent.col(vEachDim);  
    }
  }
  res.col(1) = -0.5 * sum(fitDerivError % KinvfitDerivError).t() / priorTemperature(0);
  res.col(2) = -0.5 * sum(xlatent % CinvX).t() / priorTemperature(1);
  
  // cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat eachDimensionC2(n*pdimension+theta.size(), pdimension, fill::zeros);
  for( int vEachDim = 0; vEachDim < pdimension; vEachDim++){
    if(useBand){
      bmatvecmultT(CovAllDimensions[vEachDim].mphiBand.memptr(), 
                   KinvfitDerivError.colptr(vEachDim), 
                   &(CovAllDimensions[vEachDim].bandsize), 
                   &n, 
                   eachDimensionC2.colptr(vEachDim) + n*vEachDim);
      // negate
      eachDimensionC2.col(vEachDim).subvec(n*vEachDim, n*vEachDim+n-1) =
        -eachDimensionC2.col(vEachDim).subvec(n*vEachDim, n*vEachDim+n-1);
    }else{
      eachDimensionC2.col(vEachDim).subvec(n*vEachDim, n*vEachDim+n-1) =
        -(CovAllDimensions[vEachDim].mphi.t() * KinvfitDerivError.col(vEachDim));  
    }
    
    eachDimensionC2.col(vEachDim).subvec(0, n*pdimension-1) +=
      vectorise(fderivDx.slice(vEachDim).each_col() % KinvfitDerivError.col(vEachDim));
    
    eachDimensionC2.col(vEachDim).subvec(n*pdimension, eachDimensionC2.n_rows-1) =
      fderivDtheta.slice(vEachDim).t() * KinvfitDerivError.col(vEachDim);
  }
  
  if(sigmaIsScaler){
    ret.gradient.set_size(n*pdimension + theta.size() + 1);
    ret.gradient(ret.gradient.size() - 1) = sum(sigmaGradient);
  }else{
    ret.gradient.set_size(n*pdimension + theta.size() + sigma.size());
    ret.gradient.subvec(ret.gradient.size() - sigma.size(), ret.gradient.size() - 1) = sigmaGradient;
  }
  
  ret.gradient.subvec(0, eachDimensionC2.n_rows-1) = -sum(eachDimensionC2, 1) / priorTemperature(0);
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(CinvX) / priorTemperature(1);
  fitLevelError.each_row() /= sigmaSq.t();
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(fitLevelError);
  
  return ret;
}


//' R wrapper for xthetasigmallik
//' @export
// [[Rcpp::export]]
Rcpp::List xthetasigmallikRcpp( const arma::mat & xlatent, 
                                const arma::vec & theta, 
                                const arma::vec & sigma, 
                                const arma::mat & yobs, 
                                const Rcpp::List & covAllDimInput,
                                const Rcpp::NumericVector & priorTemperatureInput = 1.0,
                                const bool useBand = false,
                                const bool useMean = false,
                                const std::string modelName = "FN"){
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'HIV'");
  }
  
  vector<gpcov> covAllDimensions(yobs.n_cols);
  for(unsigned j = 0; j < yobs.n_cols; j++){
    covAllDimensions[j] = cov_r2cpp(covAllDimInput[j]);
  }
  
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  lp ret = xthetasigmallik(xlatent, 
                           theta, 
                           sigma, 
                           yobs, 
                           covAllDimensions,
                           model,
                           priorTemperature,
                           useBand,
                           useMean);
  
  return Rcpp::List::create(Rcpp::Named("value")=ret.value,
                            Rcpp::Named("grad")=ret.gradient);
}

//' sample from GP ODE for latent x, theta, and sigma
//' @export
// [[Rcpp::export]]
Rcpp::List xthetasigmaSample( const arma::mat & yobs, 
                              const Rcpp::List & covList, 
                              const arma::vec & sigmaInit,
                              const arma::vec & xthetaInit, 
                              const arma::vec & step,
                              const int nsteps = 1, 
                              const bool traj = false, 
                              const std::string loglikflag = "usual",
                              const Rcpp::NumericVector & priorTemperatureInput = 1.0, 
                              const std::string modelName = "FN"){
  
  const arma::vec priorTemperature = Rcpp::as<arma::vec>(priorTemperatureInput);
  
  vec sigma( yobs.n_cols);
  if(sigmaInit.size() == 1){
    sigma.fill( as_scalar( sigmaInit));
  }else if(sigmaInit.size() == yobs.n_cols){
    sigma = sigmaInit;
  }else{
    throw std::runtime_error("sigmaInit size not right");
  }
  
  vector<gpcov> covAllDimensions(covList.size());
  for(unsigned int i = 0; i < covList.size(); i++){
    covAllDimensions[i] = cov_r2cpp(covList[i]);
  }
  
  OdeSystem model;
  if(modelName == "FN"){
    model = OdeSystem(fnmodelODE, fnmodelDx, fnmodelDtheta, zeros(3), ones(3)*datum::inf);
  }else if(modelName == "Hes1"){
    model = OdeSystem(hes1modelODE, hes1modelDx, hes1modelDtheta, zeros(7), ones(7)*datum::inf); 
  }else if(modelName == "HIV"){
    model = OdeSystem(HIVmodelODE, HIVmodelDx, HIVmodelDtheta, {-datum::inf, 0,0,0,0,0, -datum::inf,-datum::inf,-datum::inf}, ones(9)*datum::inf);   
  }else{
    throw std::runtime_error("modelName must be one of 'FN', 'Hes1', 'HIV'");
  }
  
  bool useBand = false;
  if(loglikflag == "band" || loglikflag == "withmeanBand"){
    useBand = true;
  }
  bool useMean = false;
  if(loglikflag == "withmean" || loglikflag == "withmeanBand"){
    useMean = true;
  }
  
  std::function<lp(vec)> tgt = [&](const vec & xthetasigma) -> lp{
      const mat & xlatent = mat(const_cast<double*>( xthetasigma.memptr()), yobs.n_rows, yobs.n_cols, false, false);
      const vec & theta = vec(const_cast<double*>( xthetasigma.memptr() + yobs.size()), xthetaInit.size() - yobs.size(), false, false);
      const vec & sigma = vec(const_cast<double*>( xthetasigma.memptr() + yobs.size() + theta.size()), sigmaInit.size(), false, false);
      return xthetasigmallik( xlatent, 
                              theta, 
                              sigma, 
                              yobs, 
                              covAllDimensions,
                              model,
                              priorTemperature,
                              useBand,
                              useMean);
    };
  
  vec lb(xthetaInit.size() + sigmaInit.size());
  lb.subvec(0, yobs.size()-1).fill(-datum::inf);
  lb.subvec(yobs.size(), xthetaInit.size()-1) = model.thetaLowerBound;
  lb.subvec(xthetaInit.size(), xthetaInit.size() + sigmaInit.size() - 1).fill(1e-3);
  
  vec initial = join_vert(xthetaInit, sigmaInit);
  hmcstate post = basic_hmcC(tgt, initial, step, lb, {datum::inf}, nsteps, traj);
  
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
