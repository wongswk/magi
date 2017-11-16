#include <cmath>
#include <random>
#include <vector>

#include "classDefinition.h"

using namespace std;

gpcov maternCov( const arma::vec &, const arma::mat &, int);
gpcov rbfCov( const arma::vec &, const arma::mat &, int);
gpcov compact1Cov( const arma::vec &, const arma::mat &, int);
lp phisigllik( const arma::vec &, const arma::mat &, const arma::mat &, string kernel = "matern");
lp xthetallik( const arma::vec & xtheta,
               const std::vector<gpcov> & CovAllDimensions,
               const double & sigma,
               const arma::mat & yobs,
               const OdeSystem & fOdeModel,
               const bool useBand = false);

lp xthetallikWithmuBand( const arma::vec & xtheta, 
                         const std::vector<gpcov> & CovAllDimensions,
                         const double & sigma, 
                         const arma::mat & yobs,
                         const OdeSystem & fOdeModel,
                         const bool useBand = true);
