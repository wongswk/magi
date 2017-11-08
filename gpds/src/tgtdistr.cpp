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
               const gpcov & CovV, 
               const gpcov & CovR, 
               const double & sigma, 
               const mat & yobs, 
               const std::function<mat (vec, mat)> & fODE) {
  int n = (xtheta.size() - 3)/2;
  const vec & theta = xtheta.subvec(xtheta.size() - 3, xtheta.size() - 1);
  lp ret;
  
  if (min(theta) < 0) {
    ret.value = -1e+9;
    ret.gradient = zeros<vec>(2*n);
    ret.gradient.subvec(xtheta.size() - 3, xtheta.size() - 1).fill(1e9);
    return ret;
  }
  
  const vec & Vsm = xtheta.subvec(0, n - 1);
  const vec & Rsm = xtheta.subvec(n, 2*n - 1);
  
  const mat & fderiv = fODE(theta, join_horiz(Vsm, Rsm));
  mat res(2,3);
  
  // V 
  vec frV = (fderiv.col(0) - CovV.mphi * Vsm); // n^2 operation
  vec VsmCTrans = CovV.CeigenVec.t() * Vsm;
  // vec frV = fderiv.col(0) - CovV.mphiLeftHalf * (VsmCTrans % CovV.Ceigen1over);
  vec frVKTrans = CovV.KeigenVec.t() * frV;
  vec fitLevelErrorV = Vsm - yobs.col(0);
  fitLevelErrorV(find_nonfinite(fitLevelErrorV)).fill(0.0);
  res(0,0) = -0.5 * sum(square( fitLevelErrorV )) / pow(sigma,2);
  res(0,1) = -0.5 * sum( square(frVKTrans) % CovV.Keigen1over);
  res(0,2) = -0.5 * sum( square(VsmCTrans) % CovV.Ceigen1over);
  
  // R
  vec frR = (fderiv.col(1) - CovR.mphi * Rsm); // n^2 operation
  vec RsmCTrans = CovR.CeigenVec.t() * Rsm;
  // vec frR = fderiv.col(1) - CovR.mphiLeftHalf * (RsmCTrans % CovR.Ceigen1over);
  vec frRKTrans = CovR.KeigenVec.t() * frR;
  vec fitLevelErrorR = Rsm - yobs.col(1);
  fitLevelErrorR(find_nonfinite(fitLevelErrorR)).fill(0.0);
  
  res(1,0) = -0.5 * sum(square( fitLevelErrorR )) / pow(sigma,2);
  res(1,1) = -0.5 * sum( square(frRKTrans) % CovR.Keigen1over);
  res(1,2) = -0.5 * sum( square(RsmCTrans) % CovR.Ceigen1over);
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = -CovV.mphi;
  Vtemp.diag() += theta(2)*(1 - square(Vsm));
  
  vec KinvFrV = (CovV.KeigenVec * (frVKTrans % CovV.Keigen1over));
  vec abcTemp = zeros<vec>(3);
  abcTemp(2) = sum(KinvFrV % fderiv.col(0)) / theta(2);
  vec VC2 =  2.0 * join_vert(join_vert( Vtemp.t()*KinvFrV, // n^2 operation
                                        theta(2) * KinvFrV ),
                                        abcTemp );
  
  
  // R contrib
  mat Rtemp = -CovR.mphi;
  Rtemp.diag() -= theta(1)/theta(2);
  
  vec KinvFrR = (CovR.KeigenVec * (frRKTrans % CovR.Keigen1over));
  abcTemp.fill(0);
  abcTemp(0) = sum(KinvFrR) / theta(2);
  abcTemp(1) = -sum(Rsm % KinvFrR) / theta(2);
  abcTemp(2) = -sum(fderiv.col(1) % KinvFrR) / theta(2);
  vec RC2 = 2.0 * join_vert(join_vert( -KinvFrR / theta(2),
                                       Rtemp.t() * KinvFrR), // n^2 operation
                                       abcTemp );
  
  vec C3 = join_vert(join_vert( 2.0 * CovV.CeigenVec * (VsmCTrans % CovV.Ceigen1over),  
                                2.0 * CovR.CeigenVec * (RsmCTrans % CovR.Ceigen1over) ), 
                                zeros<vec>(theta.size()));
  vec C1 = join_vert(join_vert( 2.0 * fitLevelErrorV / pow(sigma,2) ,  
                                2.0 * fitLevelErrorR / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}

//' log likelihood for latent states and ODE theta conditional on phi sigma
//'   latent states have a mean mu
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetallik_withmu( const vec & xtheta, 
                      const gpcov & CovV, 
                      const gpcov & CovR, 
                      const double & sigma, 
                      const mat & yobs, 
                      const OdeSystem & fOdeModel) {
  int n = (xtheta.size() - 3)/2;
  const vec & theta = xtheta.subvec(xtheta.size() - 3, xtheta.size() - 1);
  lp ret;
  
  const vec & Vsm = xtheta.subvec(0, n - 1);
  const vec & Rsm = xtheta.subvec(n, 2*n - 1);
  const vec & Vsmminusmu = Vsm - CovV.mu;
  const vec & Rsmminusmu = Rsm - CovR.mu;
  
  
  const mat & fderiv = fOdeModel.fOde(theta, join_horiz(Vsm, Rsm));
  const cube & fderivDx = fOdeModel.fOdeDx(theta, join_horiz(Vsm, Rsm));
  const cube & fderivDtheta = fOdeModel.fOdeDtheta(theta, join_horiz(Vsm, Rsm));
  
  mat res(2,3);
  
  // V 
  vec frVminusdotmu = (fderiv.col(0) - CovV.dotmu - CovV.mphi * Vsmminusmu);
  vec fitLevelErrorV = Vsm - yobs.col(0);
  fitLevelErrorV(find_nonfinite(fitLevelErrorV)).fill(0.0);
  res(0,0) = -0.5 * sum(square( fitLevelErrorV )) / pow(sigma,2);
  res(0,1) = -0.5 * as_scalar( frVminusdotmu.t() * CovV.Kinv * frVminusdotmu);
  res(0,2) = -0.5 * as_scalar( Vsmminusmu.t() * CovV.Cinv * Vsmminusmu);
  
  // R
  vec frRminusdotmu = (fderiv.col(1) - CovR.dotmu - CovR.mphi * Rsmminusmu);
  vec fitLevelErrorR = Rsm - yobs.col(1);
  fitLevelErrorR(find_nonfinite(fitLevelErrorR)).fill(0.0);
  
  res(1,0) = -0.5 * sum(square( fitLevelErrorR )) / pow(sigma,2);
  res(1,1) = -0.5 * as_scalar( frRminusdotmu.t() * CovR.Kinv * frRminusdotmu);
  res(1,2) = -0.5 * as_scalar( Rsmminusmu.t() * CovR.Cinv * Rsmminusmu);
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = -CovV.mphi;
  Vtemp.diag() += fderivDx.slice(0).col(0);
  
  vec KinvFrVminusdotmu = (CovV.Kinv * frVminusdotmu);
  
  vec VC2 =  2.0 * join_vert(join_vert( Vtemp.t() * KinvFrVminusdotmu, // n^2 operation
                                        fderivDx.slice(0).col(1) % KinvFrVminusdotmu ),
                                        fderivDtheta.slice(0).t() * KinvFrVminusdotmu );
  
  
  // R contrib
  mat Rtemp = -CovR.mphi;
  Rtemp.diag() -= fderivDx.slice(1).col(1);
  
  vec KinvFrRminusdotmu = (CovR.Kinv * frRminusdotmu);
  vec RC2 = 2.0 * join_vert(join_vert( fderivDx.slice(1).col(0) % KinvFrRminusdotmu,
                                       Rtemp.t() * KinvFrRminusdotmu), // n^2 operation
                                       fderivDtheta.slice(1).t() * KinvFrRminusdotmu );
  // 
  // vec C3 = join_vert(join_vert( 2.0 * CovV.CeigenVec * (VsmCTrans % CovV.Ceigen1over),  
  //                               2.0 * CovR.CeigenVec * (RsmCTrans % CovR.Ceigen1over) ), 
  //                               zeros<vec>(theta.size()));
  vec C3 = join_vert(join_vert( 2.0 * CovV.Cinv * Vsmminusmu,  
                                2.0 * CovR.Cinv * Rsmminusmu ), 
                                zeros<vec>(theta.size()));  
  vec C1 = join_vert(join_vert( 2.0 * fitLevelErrorV / pow(sigma,2) ,  
                                2.0 * fitLevelErrorR / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}

//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' the scale is in fact taken out and it is the legacy version of xthetallik
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetallik_rescaled( const vec & xtheta, 
                        const gpcov & CovV, 
                        const gpcov & CovR, 
                        const double & sigma, 
                        const mat & yobs, 
                        const std::function<mat (vec, mat)> & fODE) {
  vec xlatent = xtheta.subvec(0, xtheta.size() - 4);
  vec theta = xtheta.subvec(xtheta.size() - 3, xtheta.size() - 1);
  lp ret;
  
  if (min(theta) < 0) {
    ret.value = -1e+9;
    ret.gradient = zeros<vec>(xlatent.size());
    ret.gradient.subvec(xtheta.size() - 3, xtheta.size() - 1).fill(1e9);
    return ret;
  }
  
  vec Vsm = xlatent.subvec(0, xlatent.size()/2 - 1);
  vec Rsm = xlatent.subvec(xlatent.size()/2, xlatent.size() - 1);
  int n = xlatent.size()/2;
  int nobs = 0;
  
  mat fderiv = fODE(theta, join_horiz(Vsm, Rsm));
  mat res(2,3);
  
  // V 
  vec frV = (fderiv.col(0) - CovV.mphi * Vsm);
  //res(0,0) = -0.5 * sum(square( Vsm - yobs.col(0) )) / pow(sigma,2);
  res(0,0) = 0.0;
  for (int i=0; i < n; i++) {
    if (!std::isnan(yobs(i,0))) {
      res(0,0) += pow(Vsm[i] - yobs(i,0),2);
	  nobs++;
	}
  }
  res(0,0) = -0.5 * res(0,0) / pow(sigma,2);
  
  res(0,1) = -0.5 * as_scalar( frV.t() * CovV.Kinv * frV) * (double)nobs/(double)n;
  res(0,2) = -0.5 * as_scalar( Vsm.t() * CovV.Cinv * Vsm) * (double)nobs/(double)n;
  // R
  vec frR = (fderiv.col(1) - CovR.mphi * Rsm);
  res(1,0) = 0.0;
  for (int i=0; i < n; i++) {
    if (!std::isnan(yobs(i,1)))
      res(1,0) += pow(Rsm[i] - yobs(i,1),2);
  }
  res(1,0) = -0.5 * res(1,0) / pow(sigma,2);
  
  //res(1,0) = -0.5 * sum(square( Rsm - yobs.col(1) )) / pow(sigma,2);
  res(1,1) = -0.5 * as_scalar( frR.t() * CovR.Kinv * frR) * (double)nobs/(double)n;
  res(1,2) = -0.5 * as_scalar( Rsm.t() * CovR.Cinv * Rsm) * (double)nobs/(double)n;
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = eye<mat>(n, n);
  Vtemp.diag() = theta(2)*(1 - square(Vsm));
  Vtemp = Vtemp - CovV.mphi;
  mat Rtemp = eye<mat>(n, n)*theta(2);
  vec aTemp = zeros<vec>(n); 
  vec bTemp = zeros<vec>(n); 
  vec cTemp = fderiv.col(0) / theta(2);
  mat VC2 = join_horiz(join_horiz(join_horiz(join_horiz(Vtemp,Rtemp),aTemp),bTemp),cTemp);
  VC2 = 2.0 * VC2.t() * CovV.Kinv * frV * (double)nobs/(double)n;
  
  // cout << "VC2 = \n" << VC2 << endl;
  
  // R contrib
  Vtemp = eye<mat>(n, n) * -1.0/theta(2);
  Rtemp = eye<mat>(n, n) * -1.0*theta(1)/theta(2) - CovR.mphi;
  aTemp = ones<vec>(n) / theta(2);
  bTemp = -Rsm/theta(2);
  cTemp = -fderiv.col(1) / theta(2);
  mat RC2 = join_horiz(join_horiz(join_horiz(join_horiz(Vtemp,Rtemp),aTemp),bTemp),cTemp);
  RC2 = 2.0 * RC2.t() * CovR.Kinv * frR * (double)nobs/(double)n;
  
  // cout << "RC2 = \n" << RC2 << endl;
  
  vec C3 = join_vert(join_vert( 2.0 * CovV.Cinv * Vsm,  
                                2.0 * CovR.Cinv * Rsm ), 
                                zeros<vec>(theta.size()));
  C3 = C3 * (double)nobs/(double)n;
  vec C1 = join_vert(join_vert( 2.0 * (Vsm - yobs.col(0)) / pow(sigma,2) ,  
                                2.0 * (Rsm - yobs.col(1)) / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  for (unsigned int i=0; i < C1.size(); i++) {
    if (std::isnan(C1(i)))
      C1(i) = 0.0;
  }
      
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}

//' approximate log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' band matrix approximation
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
//' @noRd
//' FIXME xtheta currently passed by value for Fortran code
lp xthetallikBandApprox( const vec & xtheta, 
                         const gpcov & CovV, 
                         const gpcov & CovR, 
                         const double & sigma, 
                         const mat & yobs,
                         const std::function<mat (vec, mat)> & fODE) {
  int n = (xtheta.size() - 3)/2;
  lp ret;
  ret.gradient.set_size(xtheta.size());
  
  const double *xthetaPtr = xtheta.memptr();
  const double *VmphiPtr = CovV.mphiBand.memptr();
  const double *VKinvPtr = CovV.KinvBand.memptr();
  const double *VCinvPtr = CovV.CinvBand.memptr();
  const double *RmphiPtr = CovR.mphiBand.memptr(); 
  const double *RKinvPtr = CovR.KinvBand.memptr(); 
  const double *RCinvPtr = CovR.CinvBand.memptr();
  const double *sigmaPtr = &sigma; 
  const double *yobsPtr = yobs.memptr();
  double *retPtr = &ret.value;
  double *retgradPtr = ret.gradient.memptr();
  
  xthetallikBandC( xthetaPtr, VmphiPtr, VKinvPtr, VCinvPtr,
                   RmphiPtr, RKinvPtr, RCinvPtr, &CovV.bandsize, &n,
                   sigmaPtr, yobsPtr, retPtr, retgradPtr, fODE);
    
  return ret;
}
