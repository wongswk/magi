#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include "hmc.h"

using namespace std;
using namespace arma;

struct gpcov {
  mat C, Cinv, Kphi, Kinv, dCdphi1, dCdphi2;
};

//' matern variance covariance matrix with derivatives
//' 
//' @param phi         the parameter of (sigma_c_sq, alpha)
//' @param dist        distance matrix
//' @param complexity  how much derivative information should be calculated
gpcov maternCov(vec phi, mat dist, int complexity = 0){
  gpcov out;
  mat dist2 = square(dist);
  out.C = phi(0) * (1.0 + ((sqrt(5.0)*dist)/phi(1)) + 
    ((5.0*dist2)/(3.0*pow(phi(1),2)))) % exp((-sqrt(5.0)*dist)/phi(1));
  cout << out.C << endl;
  if (complexity == 0) return out;
  
  out.dCdphi1 = out.C/phi(0);
  out.dCdphi2 = phi(0) * ( - ((sqrt(5.0)*dist)/pow(phi(1),2)) - 
    ((10.0*dist2)/(3.0*pow(phi(1),3)))) * exp((-sqrt(5.0)*dist)/phi(1)) + 
    out.C * (sqrt(5.0)*dist)/pow(phi(1),2);
  cout << out.dCdphi2 << endl;
  if (complexity == 1) return out;
  // work from here continue for gp derivative
  return out;
}

//' log likelihood for Gaussian Process marginal likelihood with Matern kernel
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp phisigllik (vec phisig, mat yobs, mat dist){
  int n = yobs.n_rows;
  double sigma = phisig(4);
  vec res(2);
  
  // likelihood value part
  
  // V 
  gpcov CovV = maternCov(phisig.subvec(0,1), dist, 1);
  mat Kv = CovV.C+ eye<mat>(n,n)*pow(sigma, 2);
  mat Kvl = chol(Kv, "lower");
  mat Kvlinv = inv(Kvl);
  vec veta = Kvlinv * yobs.col(0);
  res(0) = -n/2.0*log(2.0*datum::pi) - sum(log(Kvl.diag())) - 0.5*sum(square(veta));
  
  // R
  gpcov CovR = maternCov(phisig.subvec(0,1), dist, 1);
  mat Kr = CovR.C+ eye<mat>(n,n)*pow(sigma, 2);
  mat Krl = chol(Kr, "lower");
  mat Krlinv = inv(Krl);
  vec reta = Krlinv * yobs.col(0);
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
  mat alphaR = Krlinv.t() * veta;
  mat facRtemp = alphaR * alphaR.t() - Krinv;
  double dRdsig = sigma * sum(facRtemp.diag());
  double dRdphi1 = accu(facRtemp % CovR.dCdphi1)/2.0;
  double dRdphi2 = accu(facRtemp % CovR.dCdphi2)/2.0;
  
  ret.gradient(0) = dVdphi1;
  ret.gradient(1) = dVdphi2;
  ret.gradient(2) = dRdphi1;
  ret.gradient(3) = dRdphi2;
  ret.gradient(4) = dVdsig+dRdsig;
  
  return ret;
}
  
  