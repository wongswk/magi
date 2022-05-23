// [[Rcpp::depends(BH)]]
#define NDEBUG
#define BOOST_DISABLE_ASSERTS

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
  // std::cout << out.C << endl;
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
    if(x < 1e-10){
        return INFINITY;
    }
  return boost::math::cyl_bessel_k(nu, x);
}

//' matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov generalMaternCov( const vec & phi,
                        const mat & distSigned,
                        int complexity){
  double noiseInjection = 1e-7;
  double df = 3.01;
  gpcov out;
  out.C.set_size(distSigned.n_rows, distSigned.n_cols);
  mat x4bessel = sqrt(2.0 * df) * abs(distSigned) / phi(1);
  
  mat bessel_df(x4bessel.n_rows, x4bessel.n_cols);
  for(unsigned int j = 0; j < bessel_df.n_cols; j++){
    for(unsigned int i = 0; i < j; i ++){
      bessel_df(i, j) = modifiedBessel2ndKind(df, x4bessel(i, j));
    }
  }
  bessel_df.diag().fill(datum::inf);
  bessel_df = symmatu(bessel_df);

  mat bessel_dfMinus1(x4bessel.n_rows, x4bessel.n_cols);
  for(unsigned int j = 0; j < bessel_df.n_cols; j++){
    for(unsigned int i = 0; i < j; i ++){
      bessel_dfMinus1(i, j) = modifiedBessel2ndKind(df-1, x4bessel(i, j));
    }
  }
  bessel_dfMinus1.diag().fill(datum::inf);
  bessel_dfMinus1 = symmatu(bessel_dfMinus1);
  
  mat bessel_dfPlus1 = bessel_dfMinus1 + 2 * df / x4bessel % bessel_df;
  bessel_dfPlus1.diag().fill(datum::inf);
  
  mat bessel_dfPlus2 = bessel_df + 2 * (df + 1) / x4bessel % bessel_dfPlus1;
  bessel_dfPlus2.diag().fill(datum::inf);

  mat bessel_dfMinus2 = bessel_df - 2 * (df - 1) / x4bessel % bessel_dfMinus1;
  bessel_dfMinus2.diag().fill(datum::inf);
  
  mat bessel_dfMinus3 = bessel_dfMinus1 - 2 * (df - 2) / x4bessel % bessel_dfMinus2;
  bessel_dfMinus3.diag().fill(datum::inf);
  
  mat bessel_dfPlus3 = bessel_dfPlus1 + 2 * (df + 2) / x4bessel % bessel_dfPlus2;
  bessel_dfPlus3.diag().fill(datum::inf);
 
  mat Cpart1 = phi(0) * pow(2.0, 1-df) * exp(-lgamma(df)) * pow( x4bessel, df);
  out.C = Cpart1 % bessel_df;
  out.C.replace(datum::nan, phi(0));
  out.C.diag() += noiseInjection;  // stabilizer

  out.mu = arma::zeros(out.C.n_rows);
  out.dotmu = arma::zeros(out.C.n_rows);
  
  if (complexity == 0) {
    return out;
  }
  
  mat dCdx4bessel = Cpart1 % (df / x4bessel % bessel_df  - 0.5 * (bessel_dfMinus1 + bessel_dfPlus1));
  dCdx4bessel.replace(datum::nan, 0);

  out.dCdphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdphiCube.slice(0) = out.C/phi(0);
  out.dCdphiCube.slice(1) = dCdx4bessel % (-sqrt(2.0 * df) / pow(phi(1), 2) * abs(distSigned));
  out.dCdphiCube.slice(1).replace(datum::nan, 0);
  
  if (complexity == 1) {
    return out;
  }
  
  // out.Cprime
  out.Cprime = dCdx4bessel % (sqrt(2.0 * df) / phi(1) * sign(distSigned));
  out.Cprime.replace(datum::nan, 0);
  
  // out.Cdoubleprime;
  mat dCprimedx4bessel =  Cpart1 * sqrt(2 * df) / phi(1);
  dCprimedx4bessel %= (
    + df * (df - 1) * pow(x4bessel, -2) % bessel_df
    - df / x4bessel % (bessel_dfMinus1 + bessel_dfPlus1)
    + 0.25 * (bessel_dfMinus2 + 2*bessel_df + bessel_dfPlus2)
  );
  dCprimedx4bessel.replace(
    datum::nan,
    -phi(0) * pow(2.0, 1-df) * exp(-lgamma(df)) * sqrt(2 * df) / phi(1) * exp(lgamma(df-1)) * pow(2, df-2)
  );
  out.Cdoubleprime = -sqrt(2 * df) / phi(1) * dCprimedx4bessel;
  
  // out.dCprimedphiCube;
  out.dCprimedphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCprimedphiCube.slice(0) = out.Cprime/phi(0);
  out.dCprimedphiCube.slice(1) = dCprimedx4bessel % (-sqrt(2.0 * df) / pow(phi(1), 2) * distSigned);  // use signed dist
  out.dCprimedphiCube.slice(1) += out.Cprime / (-phi(1));
  
  // out.dCdoubleprimedphiCube;
  out.dCdoubleprimedphiCube.set_size(out.C.n_rows, out.C.n_cols, 2);
  out.dCdoubleprimedphiCube.slice(0) = out.Cdoubleprime/phi(0);
  out.dCdoubleprimedphiCube.slice(1) = (
    + df * (df - 1) * (df - 2) * pow(x4bessel, df - 3) % bessel_df
    - 1.5 * df * (df - 1) * pow(x4bessel, df - 2) % (bessel_dfMinus1 + bessel_dfPlus1)
    + 0.75 * df * pow(x4bessel, df - 1) % (bessel_dfMinus2 + 2*bessel_df + bessel_dfPlus2)
    - 0.125 * pow(x4bessel, df) % (bessel_dfMinus3 + 3*bessel_dfMinus1 + 3*bessel_dfPlus1 + bessel_dfPlus3)
  );
  out.dCdoubleprimedphiCube.slice(1) %= phi(0) * pow(2, 1-df) * exp(-lgamma(df)) * 2 * df / pow(phi(1), 3) * x4bessel;
  out.dCdoubleprimedphiCube.slice(1) += out.Cdoubleprime * -2 / phi(1);
  const arma::uvec idx0 = arma::find(x4bessel < 1e-10);
  out.dCdoubleprimedphiCube.slice(1).elem(idx0) = out.Cdoubleprime.elem(idx0) * -2 / phi(1);
  
  // out.Cinv
  inv_sympd(out.Cinv, out.C);
  
  // out.mphi
  out.mphi = out.Cprime * out.Cinv;
  
  // out.Kinv
  out.Kphi = out.Cdoubleprime - out.mphi * out.Cprime.t();
  out.Kphi.diag() += noiseInjection;
  inv_sympd(out.Kinv, out.Kphi);
  
  // block matrix
  // TODO: for performance, I can define big matrix/cube container, and then
  // define subview<double>
  out.Sigma = join_vert(
    join_horiz(out.C, out.Cprime.t()),
    join_horiz(out.Cprime, out.Cdoubleprime)
  );
  out.dSigmadphiCube.set_size(out.Sigma.n_rows, out.Sigma.n_cols, 2);
  for(unsigned int sliceIt = 0; sliceIt < 2; sliceIt++){
    out.dSigmadphiCube.slice(sliceIt) = join_vert(
      join_horiz(out.dCdphiCube.slice(sliceIt), out.dCprimedphiCube.slice(sliceIt).t()),
      join_horiz(out.dCprimedphiCube.slice(sliceIt), out.dCdoubleprimedphiCube.slice(sliceIt))
    );
  }

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
  // std::cout << out.C << endl;
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
  // std::cout << out.C << endl;
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
  unsigned int obsDimension = yobs.n_cols;
  int phiDimension = (phisig.size() - 1) / obsDimension;
  double sigma = phisig(phisig.size() - 1);
  const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()), 
                              phiDimension, obsDimension, true, false);
  
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
    throw std::runtime_error("kernel is not specified correctly");
  }
  
  lp ret;  
  ret.gradient = zeros( phisig.size());
  ret.value = 0;
  
  // V 
  vec eigval(n);
  mat eigvec(n, n);
  
  for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
    gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
    covThisDim.C.diag() += pow(sigma, 2);
    
    eig_sym( eigval, eigvec, covThisDim.C );
    vec eta = eigvec.t() * yobs.col(pDimEach);
    ret.value += -n/2.0*log(2.0*datum::pi) - sum(log(eigval))/2.0 - 0.5*sum(square(eta) / eigval);
    
    vec alpha = eigvec * (eta / eigval);
    mat facVtemp = alpha * alpha.t() - (eigvec.each_row() % (1.0 / eigval).t()) * eigvec.t();
    // mat facVtemp = alpha * alpha.t() - inv( covThisDim.C);
    double dVdsig = sigma * sum(facVtemp.diag());
    vec dVdphi(covThisDim.dCdphiCube.n_slices);
    for(unsigned int i=0; i < dVdphi.size(); i++){
      ret.gradient(pDimEach*phiDimension + i) = accu(facVtemp % covThisDim.dCdphiCube.slice(i))/2.0;
    }
    ret.gradient(ret.gradient.size()-1) += dVdsig;
  }
  // mat CmatCholLow = chol(covThisDim.C, "lower");
  // mat CmatCholLowInv = inv(trimatl(CmatCholLow));
  // vec eta = CmatCholLowInv.t() * yobs.col(pDimEach);
  // ret.value += -n/2.0*log(2.0*datum::pi) - sum(log(CmatCholLow.diag())) - 0.5*sum(square(eta));
  // 
  // vec alpha = CmatCholLowInv.t() * eta;
  // mat facVtemp = alpha * alpha.t() - CmatCholLowInv.t() * CmatCholLowInv;
  return ret;
}

//' leave one out cross validation for Gaussian Process fitting with various kernel
//' 
//' loss function is predictive log likelihood
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigloocvllik( const vec & phisig, 
                    const mat & yobs, 
                    const mat & dist, 
                    string kernel){
  int n = yobs.n_rows;
  unsigned int obsDimension = yobs.n_cols;
  int phiDimension = (phisig.size() - 1) / obsDimension;
  const double sigma = phisig(phisig.size() - 1);
  const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()), 
                              phiDimension, obsDimension, true, false);
  
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
    throw std::runtime_error("kernel is not specified correctly");
  }
  
  lp ret;  
  ret.gradient = zeros( phisig.size());
  ret.value = 0;
  
  for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
    gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
    covThisDim.C.diag() += pow(sigma, 2);
    
    mat Cinv = arma::inv_sympd(covThisDim.C);
    
    vec alpha = Cinv * yobs.col(pDimEach);
    
    vec muLooCv = yobs.col(pDimEach) - alpha / Cinv.diag();
    vec sigmasqLooCv = 1 / Cinv.diag();
    
    ret.value += -n/2.0*log(2.0*datum::pi) - 0.5*sum(log(sigmasqLooCv)) - 
      0.5*sum(square(yobs.col(pDimEach) - muLooCv) / sigmasqLooCv);
    
    cube Ztemp = Cinv * covThisDim.dCdphiCube.each_slice();
    Ztemp = join_slices(Ztemp, 2*sigma*Cinv);
    for(unsigned int j=0; j < Ztemp.n_slices; j++){
      vec ZjCinvDiag = sum(Ztemp.slice(j).t() % Cinv).t();
      double dLdphisig = sum((alpha % (Ztemp.slice(j) * alpha) - 0.5*(1+square(alpha)/Cinv.diag())%ZjCinvDiag)/Cinv.diag());
      if(j < covThisDim.dCdphiCube.n_slices){
        ret.gradient(pDimEach*phiDimension + j) = dLdphisig;
      }else{
        ret.gradient(ret.gradient.size()-1) += dLdphisig;
      }
    }
  }
  return ret;
}

//' leave one out cross validation for Gaussian Process fitting with various kernel
//' 
//' loss function is predictive log likelihood
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigloocvmse( const vec & phisig, 
                    const mat & yobs, 
                    const mat & dist, 
                    string kernel){
  unsigned int obsDimension = yobs.n_cols;
  int phiDimension = (phisig.size() - 1) / obsDimension;
  const double sigma = phisig(phisig.size() - 1);
  const mat & phiAllDim = mat(const_cast<double*>( phisig.begin()), 
                              phiDimension, obsDimension, true, false);
  
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
    throw std::runtime_error("kernel is not specified correctly");
  }
  
  lp ret;  
  ret.gradient = zeros( phisig.size());
  ret.value = 0;
  
  for(unsigned int pDimEach = 0; pDimEach < obsDimension; pDimEach++){
    gpcov covThisDim = kernelCov(phiAllDim.col(pDimEach), dist, 1);
    covThisDim.C.diag() += pow(sigma, 2);
    
    mat Cinv = arma::inv_sympd(covThisDim.C);
    
    vec alpha = Cinv * yobs.col(pDimEach);
    
    vec muLooCv = yobs.col(pDimEach) - alpha / Cinv.diag();
    vec sigmasqLooCv = 1 / Cinv.diag();
    
    ret.value += -0.5*sum(square(yobs.col(pDimEach) - muLooCv));
    
    cube Ztemp = Cinv * covThisDim.dCdphiCube.each_slice();
    Ztemp = join_slices(Ztemp, 2*sigma*Cinv);
    for(unsigned int j=0; j < Ztemp.n_slices; j++){
      vec ZjCinvDiag = sum(Ztemp.slice(j).t() % Cinv).t();
      vec dmudphisig = (Ztemp.slice(j) * alpha) / Cinv.diag() - alpha % ZjCinvDiag / square(Cinv.diag());
      if(j < covThisDim.dCdphiCube.n_slices){
        ret.gradient(pDimEach*phiDimension + j) = sum(dmudphisig % (yobs.col(pDimEach) - muLooCv));
      }else{
        ret.gradient(ret.gradient.size()-1) += sum(dmudphisig % (yobs.col(pDimEach) - muLooCv));
      }
    }
  }
  return ret;
}

//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetallik( const vec & xtheta, 
               const std::vector<gpcov> & CovAllDimensions, 
               const vec & sigma, 
               const mat & yobs, 
               const OdeSystem & fOdeModel,
               const bool useBand,
               const arma::vec & priorTemperatureInput) {
  const arma::vec & tvecFull = CovAllDimensions[0].tvecCovInput;
  int n = yobs.n_rows;
  int pdimension = yobs.n_cols;
  const mat & xlatent = mat(const_cast<double*>( xtheta.memptr()), n, pdimension, false, false);
  const vec & theta = xtheta.subvec(n*pdimension, xtheta.size() - 1);
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

  const mat & fderiv = fOdeModel.fOde(theta, xlatent, tvecFull);
  const cube & fderivDx = fOdeModel.fOdeDx(theta, xlatent, tvecFull);
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, xlatent, tvecFull);
  
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
  res.col(0) = -0.5 * sum(square( fitLevelError )).t() / square(sigma);
  res.col(0) /= priorTemperature(2);
  
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
  
  // std::cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // std::cout << "lglik = " << ret.value << endl;
  
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
  
  ret.gradient = -sum(eachDimensionC2, 1) / priorTemperature(0);
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(CinvX) / priorTemperature(1);
  rowvec sigmaSq = square(sigma.t());
  fitLevelError.each_row() /= sigmaSq;
  fitLevelError /= priorTemperature(2);
  ret.gradient.subvec(0, n*pdimension-1) -= vectorise(fitLevelError);
  
  return ret;
}


// log likelihood for latent states and ODE theta conditional on phi sigma
// with mean 
lp xthetallikWithmuBand( const vec & xtheta, 
                         const std::vector<gpcov> & CovAllDimensions,
                         const vec & sigma, 
                         const mat & yobs, 
                         const OdeSystem & fOdeModel,
                         const bool useBand,
                         const arma::vec & priorTemperature) {
  int n = yobs.n_rows;
  vec xthetaShifted = xtheta;
  mat yobsShifted = yobs;
  mat muAllDimension(yobs.n_rows, yobs.n_cols);
  mat dotmuAllDimension(yobs.n_rows, yobs.n_cols);
  
  for(unsigned int i = 0; i < yobs.n_cols; i++){
    xthetaShifted.subvec(i*n, i*n + n - 1) -= CovAllDimensions[i].mu;
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
  
  lp ret = xthetallik(xthetaShifted, CovAllDimensions, sigma, yobsShifted, 
                      fOdeModelShifted, useBand, priorTemperature); 
  return ret;
}

