#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
#include "hmc.h"
#include "dynamicalSystemModels.h"
#include "band.h"


using namespace arma;


//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetaphi1sigmallik( const mat & xlatent, 
                        const vec & theta, 
                        const vec & phi1, 
                        const vec & sigmaInput, 
                        const mat & yobs, 
                        const std::vector<gpcov> & CovAllDimensions,
                        const OdeSystem & fOdeModel,
                        const arma::vec & priorTemperatureInput,
                        const bool useBand,
                        const bool useMean) {
  const arma::vec & tvecFull = CovAllDimensions[0].tvecCovInput;
  
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
    (const vec & theta, const mat & x, const vec & tvec) -> mat{
      return fOdeModel.fOde(theta, x+muAllDimension, tvec) - dotmuAllDimension;
    };
    
    fOdeModelShifted.fOdeDx = [&muAllDimension, &dotmuAllDimension, &fOdeModel]
    (const vec & theta, const mat & x, const vec & tvec) -> cube{
      return fOdeModel.fOdeDx(theta, x+muAllDimension, tvec);
    };
    
    fOdeModelShifted.fOdeDtheta = [&muAllDimension, &dotmuAllDimension, &fOdeModel]
    (const vec & theta, const mat & x, const vec & tvec) -> cube{
      return fOdeModel.fOdeDtheta(theta, x+muAllDimension, tvec);
    };
    
    return xthetaphi1sigmallik(xlatentShifted, theta, phi1, sigmaInput, yobsShifted, CovAllDimensions, 
                               fOdeModelShifted, priorTemperatureInput, useBand, false); 
  }
  
  int n = yobs.n_rows;
  int pdimension = yobs.n_cols;
  
  lp ret;
  if (fOdeModel.checkBound(xlatent, theta, &ret)) {
    return ret;

  }

  arma::vec priorTemperature(3);
  if(priorTemperatureInput.n_rows == 1){
    priorTemperature.fill(as_scalar(priorTemperatureInput));
  }else if(priorTemperatureInput.n_rows == 2){
    priorTemperature.subvec(0, 1) = priorTemperatureInput;
    priorTemperature(2) = 1.0;
  }else if(priorTemperatureInput.n_rows == 3){
    priorTemperature = priorTemperatureInput;
  }else{
    throw std::invalid_argument("priorTemperatureInput must be scaler, 2-vector or 3-vector");
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
  
  const mat & fderiv = fOdeModel.fOde(theta, xlatent, tvecFull);
  const cube & fderivDx = fOdeModel.fOdeDx(theta, xlatent, tvecFull);
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, xlatent, tvecFull);
  
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
  res.col(0) /= priorTemperature(2);
  vec sigmaGradient = sum(square( fitLevelError )).t() / (sigmaSq % sigma) - nobs / sigma;
  sigmaGradient /= priorTemperature(2);
  
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
  res.col(1) = -0.5 * sum(fitDerivError % KinvfitDerivError).t() / phi1 / priorTemperature(0);
  res.col(2) = -0.5 * sum(xlatent % CinvX).t() / phi1 / priorTemperature(1);
  
  // cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  ret.value -= xlatent.n_rows * accu(log(phi1));
  
  // for(unsigned int pDimEach = 0; pDimEach < yobs.n_cols; pDimEach++){
  //   const gpcov & covThisDim = CovAllDimensions[pDimEach];
  //   ret.value += -xlatent.n_rows*log(2.0*datum::pi) - 0.5*covThisDim.logCdet - 0.5*covThisDim.logKdet;
  // }
  
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
    ret.gradient.set_size(n*pdimension + theta.size() + phi1.size() + 1);
    ret.gradient(ret.gradient.size() - 1) = sum(sigmaGradient);
  }else{
    ret.gradient.set_size(n*pdimension + theta.size() + phi1.size() + sigma.size());
    ret.gradient.subvec(ret.gradient.size() - sigma.size(), ret.gradient.size() - 1) = sigmaGradient;
  }
  
  ret.gradient.subvec(0, eachDimensionC2.n_rows-1) = -sum(eachDimensionC2.each_row() / phi1.t(), 1) / priorTemperature(0);
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(CinvX.each_row() / phi1.t()) / priorTemperature(1);
  fitLevelError.each_row() /= sigmaSq.t();
  fitLevelError /= priorTemperature(2);
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(fitLevelError);
  
  ret.gradient.subvec(eachDimensionC2.n_rows, eachDimensionC2.n_rows + phi1.size() - 1) = 
    (-res.col(1) / phi1) + (-res.col(2) / phi1) - xlatent.n_rows / phi1;
  
  return ret;
}
