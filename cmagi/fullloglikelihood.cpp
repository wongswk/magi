#include "classDefinition.h"
#include "fullloglikelihood.h"
#include "tgtdistr.h"
#include "hmc.h"
#include "dynamicalSystemModels.h"

using namespace arma;


//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetaphisigmallik( const mat & xlatent, 
                       const vec & theta, 
                       const mat & phi, 
                       const vec & sigmaInput, 
                       const mat & yobs, 
                       const vec & xtimes,
                       const OdeSystem & fOdeModel) {
  mat distSigned(xtimes.size(), xtimes.size());
  for(unsigned int i = 0; i < distSigned.n_cols; i++){
    distSigned.col(i) = xtimes - xtimes(i);
  }
  
  vector<gpcov> CovAllDimensions(phi.n_cols);
  for(unsigned int j = 0; j < phi.n_cols; j++){
    CovAllDimensions[j] = generalMaternCov( phi.col(j), distSigned, 3);
  }
  
  vec xtheta = join_vert(vectorise(xlatent), theta);
  
  vec sigma(phi.n_cols);
  bool sigmaIsScaler = (sigmaInput.size() == 1);
  if (sigmaIsScaler){
    sigma.fill(as_scalar(sigmaInput));
  }else if(sigmaInput.size() == phi.n_cols){
    sigma = sigmaInput;
  }else{
    throw std::runtime_error("sigmaInput dimension not right");
  }
  
  lp numerator = xthetallik(xtheta, CovAllDimensions, sigma, yobs, fOdeModel, false, ones(2));
  
  const mat & fderiv = fOdeModel.fOde(theta, xlatent);
  mat phiGradient(phi.n_rows, phi.n_cols);
  vec sigmaGradient(sigma.size());
  
  lp bigpost;
  for(unsigned int pDimEach = 0; pDimEach < yobs.n_cols; pDimEach++){
    const gpcov & covThisDim = CovAllDimensions[pDimEach];
    vec eigval(covThisDim.Sigma.n_cols);
    mat eigvec(covThisDim.Sigma.n_rows, covThisDim.Sigma.n_cols);
    
    eig_sym( eigval, eigvec, covThisDim.Sigma );
    vec eta = eigvec.t() * join_vert(xlatent.col(pDimEach), fderiv.col(pDimEach));
    
    bigpost.value += -0.5*covThisDim.Sigma.n_rows*log(2.0*datum::pi) - sum(log(eigval))/2.0 - 0.5*sum(square(eta) / eigval);
    vec fitLevelError = yobs.col(pDimEach) - xlatent.col(pDimEach);
    int nobs = uvec(find_finite(yobs.col(pDimEach))).size();
    fitLevelError(find_nonfinite(fitLevelError)).fill(0.0);
    bigpost.value += -0.5*log(2.0*datum::pi)*nobs - log(sigma(pDimEach))*nobs - 0.5 * sum(square(fitLevelError) / pow(sigma(pDimEach), 2));
    
    vec alpha = eigvec * (eta / eigval);
    mat facVtemp = alpha * alpha.t() - (eigvec.each_row() % (1.0 / eigval).t()) * eigvec.t();
    
    for(unsigned int i=0; i < covThisDim.dSigmadphiCube.n_slices; i++){
      phiGradient(i, pDimEach) = 0.5*accu(facVtemp % covThisDim.dSigmadphiCube.slice(i));
    }
    
    sigmaGradient(pDimEach) = sum(square(fitLevelError)) / pow(sigma(pDimEach), 3) - nobs / sigma(pDimEach);
  }
  
  if(sigmaIsScaler){
    bigpost.gradient.set_size(numerator.gradient.size() + phiGradient.size() + 1);
    bigpost.gradient(bigpost.gradient.size()-1) = sum(sigmaGradient);
  }else{
    bigpost.gradient.set_size(numerator.gradient.size() + phiGradient.size() + sigma.size());
    bigpost.gradient.subvec(numerator.gradient.size() + phiGradient.size(), bigpost.gradient.size()-1) = sigmaGradient;
  }
  bigpost.gradient.subvec(0, numerator.gradient.size()-1) = numerator.gradient;
  bigpost.gradient.subvec(numerator.gradient.size(), numerator.gradient.size() + phiGradient.size()-1) = vectorise(phiGradient);
  return bigpost;
}
