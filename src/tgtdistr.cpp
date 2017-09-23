#include "tgtdistr.h"
#include "hmc.h"

//' matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov maternCov( vec phi, mat dist, int complexity = 0){
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * (1.0 + ((sqrt(5.0)*dist)/phi(1)) + 
    ((5.0*dist2)/(3.0*pow(phi(1),2)))) % exp((-sqrt(5.0)*dist)/phi(1));
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphi1 = out.C/phi(0);
  out.dCdphi2 = phi(0) * ( - ((sqrt(5.0)*dist)/pow(phi(1),2)) - 
    ((10.0*dist2)/(3.0*pow(phi(1),3)))) % exp((-sqrt(5.0)*dist)/phi(1)) + 
    out.C % ((sqrt(5.0)*dist)/pow(phi(1),2));
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}

gpcov rbfCov( vec phi, mat dist, int complexity = 0){
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * exp(-dist2/(2.0*pow(phi(1), 2)));
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphi1 = out.C/phi(0);
  out.dCdphi2 = out.C % dist2 / pow(phi(1), 3);
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}

gpcov compact1Cov( vec phi, mat dist, int complexity = 0){
  int dimension = 3;
  mat zeromat = zeros<mat>(dist.n_rows, dist.n_cols);
  int p = floor((double)dimension / 2.0) + 2;
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * pow( arma::max(1 - dist / phi(1), zeromat), p+1) % ((p+1)*dist/phi(1)+1);
  out.C.diag() += 1e-7;
  // cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphi1 = out.C/phi(0);
  out.dCdphi2 = phi(0) * pow( arma::max(1 - dist / phi(1), zeromat), p) 
    % pow(dist,2)/pow(phi(1),3) * (p+1) * (p+2);
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}


//' log likelihood for Gaussian Process marginal likelihood with Matern kernel
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigllik( vec phisig, 
               mat yobs, 
               mat dist, 
               string kernel){
  int n = yobs.n_rows;
  double sigma = phisig(4);
  vec res(2);
  
  // likelihood value part
  std::function<gpcov(vec, mat, int)> kernelCov;
  if(kernel == "matern"){
    kernelCov = maternCov;
  }else if(kernel == "rbf"){
    kernelCov = rbfCov;
  }else if(kernel == "compact1"){
    kernelCov = compact1Cov;
  }else{
    throw "kernel is not specified correctly";
  }
  
  // V 
  gpcov CovV = kernelCov(phisig.subvec(0,1), dist, 1);
  mat Kv = CovV.C+ eye<mat>(n,n)*pow(sigma, 2);
  mat Kvl = chol(Kv, "lower");
  mat Kvlinv = inv(Kvl);
  vec veta = Kvlinv * yobs.col(0);
  res(0) = -n/2.0*log(2.0*datum::pi) - sum(log(Kvl.diag())) - 0.5*sum(square(veta));
  
  // R
  gpcov CovR = kernelCov(phisig.subvec(2,3), dist, 1);
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
  double dVdphi1 = accu(facVtemp % CovV.dCdphi1)/2.0;
  double dVdphi2 = accu(facVtemp % CovV.dCdphi2)/2.0;
  
  // R contrib
  mat Krinv = Krlinv.t() * Krlinv;
  mat alphaR = Krlinv.t() * reta;
  mat facRtemp = alphaR * alphaR.t() - Krinv;
  double dRdsig = sigma * sum(facRtemp.diag());
  double dRdphi1 = accu(facRtemp % CovR.dCdphi1)/2.0;
  double dRdphi2 = accu(facRtemp % CovR.dCdphi2)/2.0;
  
  ret.gradient = zeros<vec>(5);
  ret.gradient(0) = dVdphi1;
  ret.gradient(1) = dVdphi2;
  ret.gradient(2) = dRdphi1;
  ret.gradient(3) = dRdphi2;
  ret.gradient(4) = dVdsig+dRdsig;
  
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
  vec frV = (fderiv.col(0) - CovV.mphi * Vsm);
  vec fitLevelErrorV = Vsm - yobs.col(0);
  fitLevelErrorV(find_nonfinite(fitLevelErrorV)).fill(0.0);
  res(0,0) = -0.5 * sum(square( fitLevelErrorV )) / pow(sigma,2);
  res(0,1) = -0.5 * as_scalar( frV.t() * CovV.Kinv * frV);
  res(0,2) = -0.5 * as_scalar( Vsm.t() * CovV.Cinv * Vsm);
  
  // R
  vec frR = (fderiv.col(1) - CovR.mphi * Rsm);
  vec fitLevelErrorR = Rsm - yobs.col(1);
  fitLevelErrorR(find_nonfinite(fitLevelErrorR)).fill(0.0);
  
  res(1,0) = -0.5 * sum(square( fitLevelErrorR )) / pow(sigma,2);
  res(1,1) = -0.5 * as_scalar( frR.t() * CovR.Kinv * frR);
  res(1,2) = -0.5 * as_scalar( Rsm.t() * CovR.Cinv * Rsm);
  
  //cout << "lglik component = \n" << res << endl;
  
  ret.value = accu(res);
  
  // cout << "lglik = " << ret.value << endl;
  
  // gradient 
  // V contrib
  mat Vtemp = eye<mat>(n, n);
  Vtemp.diag() = theta(2)*(1 - square(Vsm));
  Vtemp -= CovV.mphi;
  mat Rtemp = eye<mat>(n, n)*theta(2);
  vec aTemp = zeros<vec>(n); 
  vec bTemp = zeros<vec>(n); 
  vec cTemp = fderiv.col(0) / theta(2);
  mat VC2 = join_horiz(join_horiz(join_horiz(join_horiz(Vtemp,Rtemp),aTemp),bTemp),cTemp);
  VC2 = 2.0 * VC2.t() * (CovV.Kinv * frV);
  
  // cout << "VC2 = \n" << VC2 << endl;
  
  // R contrib
  Vtemp = zeros<mat>(n, n);
  Vtemp.diag().fill(-1.0/theta(2));
  Rtemp = zeros<mat>(n, n);
  Rtemp.diag().fill(-1.0*theta(1)/theta(2));
  Rtemp -= CovR.mphi;
  aTemp = ones<vec>(n) / theta(2);
  bTemp = -Rsm/theta(2);
  cTemp = -fderiv.col(1) / theta(2);
  mat RC2 = join_horiz(join_horiz(join_horiz(join_horiz(Vtemp,Rtemp),aTemp),bTemp),cTemp);
  RC2 = 2.0 * RC2.t() * (CovR.Kinv * frR);
  
  // cout << "RC2 = \n" << RC2 << endl;
  
  vec C3 = join_vert(join_vert( 2.0 * CovV.Cinv * Vsm,  
                                2.0 * CovR.Cinv * Rsm ), 
                                zeros<vec>(theta.size()));
  vec C1 = join_vert(join_vert( 2.0 * fitLevelErrorV / pow(sigma,2) ,  
                                2.0 * fitLevelErrorR / pow(sigma,2) ),
                                zeros<vec>(theta.size()));
  
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}

//' log likelihood for latent states and ODE theta conditional on phi sigma
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
  for (int i=0; i < C1.size(); i++) {
    if (std::isnan(C1(i)))
      C1(i) = 0.0;
  }
      
  ret.gradient = ((VC2 + RC2)  + C3 + C1 ) * -0.5;
  
  return ret;
}



mat fnmodelODE(const vec & theta, const mat & x) {
  vec V = x.col(0);
  vec R = x.col(1);
  
  vec Vdt = theta(2) * (V - pow(V,3) / 3.0 + R);
  vec Rdt = -1.0/theta(2) * ( V - theta(0) + theta(1) * R);
  
  return join_horiz(Vdt, Rdt);
}