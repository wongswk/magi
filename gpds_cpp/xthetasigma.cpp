#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
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
                    const arma::vec & priorTemperatureInput,
                    const bool useBand,
                    const bool useMean) {
  
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
  res.col(1) = -0.5 * sum(fitDerivError % KinvfitDerivError).t() / priorTemperature(0);
  res.col(2) = -0.5 * sum(xlatent % CinvX).t() / priorTemperature(1);

  if(arma::any(res.col(2) > arma::min(arma::abs(res.col(0)), arma::abs(res.col(1))))){
      std::cout << "lglik component = \n" << res << std::endl;
      throw std::runtime_error("smoothing component positive definiteness violated, consider increase band size");
  }

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
