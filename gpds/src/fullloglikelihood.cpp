#include "fullloglikelihood.h"
#include "tgtdistr.h"
// #include "dynamicalSystemModels.h"
// #include <boost/math/special_functions/bessel.hpp>

using namespace arma;


//' log likelihood for latent states and ODE theta conditional on phi sigma
//' 
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetaphisigmallik( const mat & xlatent, 
                       const vec & theta, 
                       const mat & phi, 
                       const vec & sigma, 
                       const mat & yobs, 
                       const vec & xtimes,
                       const OdeSystem & fOdeModel) {
  mat distSigned(xtimes.size(), xtimes.size());
  for(unsigned int i = 0; i < distSigned.n_cols; i++){
    distSigned.col(i) = xtimes - xtimes(i);
  }
  cout << "mark 1\n";
  vector<gpcov> CovAllDimensions(phi.n_cols);
  for(unsigned int j = 0; j < phi.n_cols; j++){
    CovAllDimensions[j] = generalMaternCov( phi.col(j), distSigned, 3);
  }
  cout << "mark 2\n";
  vec xtheta = join_vert(vectorise(xlatent), theta);
  cout << "mark 3\n";
  lp numerator = xthetallik(xtheta, CovAllDimensions, sigma, yobs, fOdeModel, false, ones(2));
  cout << "mark 4\n";
  
  return numerator;
}
