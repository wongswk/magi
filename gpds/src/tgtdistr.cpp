#include "tgtdistr.h"
#include "band.h"
#include "dynamicalSystemModels.h"
#include <boost/math/special_functions/bessel.hpp>

using namespace arma;

//' matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov maternCov( const vec & phi, const mat & dist, int complexity = 0){
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * (1.0 + ((sqrt(5.0)*dist)/phi(1)) + 
    ((5.0*dist2)/(3.0*pow(phi(1),2)))) % exp((-sqrt(5.0)*dist)/phi(1));
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdphiCube.slice(0) = out.C/phi(0);
  out.dCdphiCube.slice(1) = phi(0) * ( - ((sqrt(5.0)*dist)/pow(phi(1),2)) - 
    ((10.0*dist2)/(3.0*pow(phi(1),3)))) % exp((-sqrt(5.0)*dist)/phi(1)) + 
    out.C % ((sqrt(5.0)*dist)/pow(phi(1),2));
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}

double modifiedBessel2ndKind (const double & nu, const double & x){
  return boost::math::cyl_bessel_k(nu, x);
}

//' matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov generalMaternCov( const vec & phi, const mat & dist, int complexity = 0){
  double df = 2.01;
  gpcov out;
  mat dist2 = square(dist);
  out.C.set_size(dist.n_rows, dist.n_cols);
  mat x4bessel = sqrt(2.0 * df) * dist / phi(1);
  for(unsigned int i = 0; i < dist.size(); i++){
    if(abs(dist(i)) < 1e-14){
      out.C(i) = phi(0);
    }else{
      out.C(i) = phi(0) * pow(2.0, 1-df) * exp(-lgamma(df)) * pow( x4bessel(i), df) * 
        modifiedBessel2ndKind(df, x4bessel(i));  
    }
  }
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  out.dCdphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdphiCube.slice(0) = out.C/phi(0);
  out.dCdphiCube.slice(1) = out.C * df / x4bessel;
  for(unsigned int i = 0; i < dist.size(); i++){
    if(abs(dist(i)) < 1e-14){
      out.dCdphiCube.slice(1)(i) = 0;
    }else{
      out.dCdphiCube.slice(1)(i) += phi(0) * pow(2.0, 1-df) * exp(-lgamma(df)) * pow( x4bessel(i), df) * 
        -(modifiedBessel2ndKind(df-1, x4bessel(i)) + modifiedBessel2ndKind(df+1, x4bessel(i)))/2.0;  
    }
    out.dCdphiCube.slice(1)(i) *= -sqrt(2.0 * df) * dist(i) / pow(phi(1), 2);
  }
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}


//' periodic matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov periodicMaternCov( const vec & phi, const mat & dist, int complexity = 0){
  mat newdist = abs(sin(dist * datum::pi / phi(2))) * 2.0;
  gpcov out = maternCov( phi.subvec(0,1), newdist, complexity);
  out.dCdphiCube.resize(out.dCdphiCube.n_rows, out.dCdphiCube.n_cols, 3);
  out.dCdphiCube.slice(2) = out.C % sign(sin(dist*datum::pi/phi(2))) 
    % (cos(dist*datum::pi/phi(2))*2) % (dist*datum::pi * -1/pow(phi(2),2));
  return out;
}
gpcov rbfCov( const vec & phi, const mat & dist, int complexity = 0){
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * exp(-dist2/(2.0*pow(phi(1), 2)));
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdphiCube.slice(0) = out.C/phi(0);
  out.dCdphiCube.slice(1) = out.C % dist2 / pow(phi(1), 3);
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}

gpcov compact1Cov( const vec & phi, const mat & dist, int complexity = 0){
  int dimension = 3;
  mat zeromat = zeros<mat>(dist.n_rows, dist.n_cols);
  int p = floor((double)dimension / 2.0) + 2;
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * pow( arma::max(1 - dist / phi(1), zeromat), p+1) % ((p+1)*dist/phi(1)+1);
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdphiCube.slice(0) = out.C/phi(0);
  out.dCdphiCube.slice(1) = phi(0) * pow( arma::max(1 - dist / phi(1), zeromat), p) 
    % pow(dist,2)/pow(phi(1),3) * (p+1) * (p+2);
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}


//' log likelihood for Gaussian Process marginal likelihood with Matern kernel
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigllik( const vec & phisig, 
               const mat & yobs, 
               const mat & dist, 
               string kernel){
  int n = yobs.n_rows;
  int p = phisig.size();
  double sigma = phisig(p-1);
  vec res(2);
  
  // likelihood value part
  std::function<gpcov(vec, mat, int)> kernelCov;
  if(kernel == "matern"){
    kernelCov = maternCov;
  }else if(kernel == "rbf"){
    kernelCov = rbfCov;
  }else if(kernel == "compact1"){
    kernelCov = compact1Cov;
  }else if(kernel == "periodicMatern"){
    kernelCov = periodicMaternCov;
  }else if(kernel == "generalMatern"){
    kernelCov = generalMaternCov;
  }else{
    throw "kernel is not specified correctly";
  }
  
  // V 
  gpcov CovV = kernelCov(phisig.subvec(0,(p-1)/2-1), dist, 1);
  mat Kv = CovV.C+ eye<mat>(n,n)*pow(sigma, 2);
  mat Kvl = chol(Kv, "lower");
  mat Kvlinv = inv(Kvl);
  vec veta = Kvlinv * yobs.col(0);
  res(0) = -n/2.0*log(2.0*datum::pi) - sum(log(Kvl.diag())) - 0.5*sum(square(veta));
  
  // R
  gpcov CovR = kernelCov(phisig.subvec((p-1)/2,p-2), dist, 1);
  mat Kr = CovR.C+ eye<mat>(n,n)*pow(sigma, 2);
  mat Krl = chol(Kr, "lower");
  mat Krlinv = inv(Krl);
  vec reta = Krlinv * yobs.col(1);
  res(1) = -n/2.0*log(2.0*datum::pi) - sum(log(Krl.diag())) - 0.5*sum(square(reta));
  
  lp ret;  
  ret.value = sum(res);
  
  // gradient part
  // V contrib
  mat Kvinv = Kvlinv.t() * Kvlinv;
  mat alphaV = Kvlinv.t() * veta;
  mat facVtemp = alphaV * alphaV.t() - Kvinv;
  double dVdsig = sigma * sum(facVtemp.diag());
  
  vec dVdphi(CovV.dCdphiCube.n_slices);
  for(unsigned int i=0; i<dVdphi.size(); i++){
    dVdphi(i) = accu(facVtemp % CovV.dCdphiCube.slice(i))/2.0;
  }
  
  // R contrib
  mat Krinv = Krlinv.t() * Krlinv;
  mat alphaR = Krlinv.t() * reta;
  mat facRtemp = alphaR * alphaR.t() - Krinv;
  double dRdsig = sigma * sum(facRtemp.diag());
  
  vec dRdphi(CovR.dCdphiCube.n_slices);
  for(unsigned int i=0; i<dRdphi.size(); i++){
    dRdphi(i) = accu(facRtemp % CovR.dCdphiCube.slice(i))/2.0;
  }
  
  ret.gradient.set_size(p);
  ret.gradient.subvec(0, dVdphi.size()-1) = dVdphi;
  ret.gradient.subvec(dVdphi.size(), p-2) = dRdphi;
  ret.gradient(p-1) = dVdsig+dRdsig;
  
  // cout << ret.value << endl << ret.gradient;
  
  return ret;
}

//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetallik( const vec & xtheta, 
               const std::vector<gpcov> & CovAllDimensions, 
               const double & sigma, 
               const mat & yobs, 
               const OdeSystem & fOdeModel,
               const bool useBand,
               const double & priorTemperature) {
  int n = yobs.n_rows;
  int pdimension = yobs.n_cols;
  const mat & xlatent = mat(const_cast<double*>( xtheta.memptr()), n, pdimension, false, false);
  const vec & theta = xtheta.subvec(n*pdimension, xtheta.size() - 1);
  lp ret;
  if (fOdeModel.checkBound(xlatent, theta, &ret)) {
    return ret;
  }
  
  const mat & fderiv = fOdeModel.fOde(theta, xlatent);
  const cube & fderivDx = fOdeModel.fOdeDx(theta, xlatent);
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, xlatent);
  
  mat res(pdimension, 3);
  
  // V 
  mat fitLevelError = xlatent - yobs;
  mat fitDerivError(n, pdimension);
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
  }
  fitLevelError(find_nonfinite(fitLevelError)).fill(0.0);
  res.col(0) = -0.5 * sum(square( fitLevelError )).t() / pow(sigma,2);
  
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
  res.col(1) = -0.5 * sum(fitDerivError % KinvfitDerivError).t() / priorTemperature;
  res.col(2) = -0.5 * sum(xlatent % CinvX).t() / priorTemperature;
  
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
  
  ret.gradient = -sum(eachDimensionC2, 1) / priorTemperature;
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(CinvX) / priorTemperature;
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(fitLevelError) / pow(sigma, 2);
  
  return ret;
}


// log likelihood for latent states and ODE theta conditional on phi sigma
// with mean 
lp xthetallikWithmuBand( const vec & xtheta, 
                         const std::vector<gpcov> & CovAllDimensions,
                         const double & sigma, 
                         const mat & yobs, 
                         const OdeSystem & fOdeModel,
                         const bool useBand,
                         const double & priorTemperature) {
  const gpcov & CovV = CovAllDimensions[0];
  const gpcov & CovR = CovAllDimensions[1];
  int n = (xtheta.size() - 3)/2;
  vec xthetaShifted = xtheta;
  xthetaShifted.subvec(0, n - 1) -= CovV.mu;
  xthetaShifted.subvec(n, 2*n - 1) -= CovR.mu;
  
  mat yobsShifted = yobs;
  yobsShifted.col(0) -= CovV.mu;
  yobsShifted.col(1) -= CovR.mu;
  
  OdeSystem fOdeModelShifted = fOdeModel;
  
  fOdeModelShifted.fOde = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> mat{
    return fOdeModel.fOde(theta, x+join_horiz(CovV.mu, CovR.mu)) - join_horiz(CovV.dotmu, CovR.dotmu);
  };
  
  fOdeModelShifted.fOdeDx = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> cube{ 
    return fOdeModel.fOdeDx(theta, x+join_horiz(CovV.mu, CovR.mu));
  };
  
  fOdeModelShifted.fOdeDtheta = [&CovV, &CovR, &fOdeModel](const vec & theta, const mat & x) -> cube{ 
    return fOdeModel.fOdeDtheta(theta, x+join_horiz(CovV.mu, CovR.mu));
  };
  
  lp ret = xthetallik(xthetaShifted, CovAllDimensions, sigma, yobsShifted, fOdeModelShifted, useBand, priorTemperature); 
  return ret;
}

