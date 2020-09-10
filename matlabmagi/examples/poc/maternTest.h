#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>
#include <armadillo>

struct gpcov {
  arma::mat C, Cprime, Cdoubleprime, Cinv, mphi, Kphi, Kinv, CeigenVec, KeigenVec, mphiLeftHalf;
  arma::mat Sigma;
  arma::cube dCdphiCube, dCprimedphiCube, dCdoubleprimedphiCube, dSigmadphiCube;
  arma::mat CinvBand, mphiBand, KinvBand;
  arma::vec Ceigen1over, Keigen1over, mu, dotmu;
  int bandsize;
};

gpcov maternCov( const arma::vec & phi, const arma::mat & dist, int complexity);
