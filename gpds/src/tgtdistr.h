#include <cmath>
#include <random>
#include <vector>

#include "classDefinition.h"

using namespace std;

gpcov maternCov( const arma::vec &, const arma::mat &, int);
gpcov generalMaternCov( const arma::vec &, const arma::mat &, int);
gpcov rbfCov( const arma::vec &, const arma::mat &, int);
gpcov compact1Cov( const arma::vec &, const arma::mat &, int);
lp phisigllik( const arma::vec &, const arma::mat &, const arma::mat &, string kernel = "matern");
lp phisigloocvllik( const arma::vec &, const arma::mat &, const arma::mat &, string kernel = "matern");
lp phisigloocvmse( const arma::vec &, const arma::mat &, const arma::mat &, string kernel = "matern");
lp xthetallik( const arma::vec & xtheta,
               const std::vector<gpcov> & CovAllDimensions,
               const arma::vec & sigma,
               const arma::mat & yobs,
               const OdeSystem & fOdeModel,
               const bool useBand = false,
               const arma::vec & priorTemperatureInput = arma::ones(2));

lp xthetallikWithmuBand( const arma::vec & xtheta, 
                         const std::vector<gpcov> & CovAllDimensions,
                         const arma::vec & sigma, 
                         const arma::mat & yobs,
                         const OdeSystem & fOdeModel,
                         const bool useBand = true,
                         const arma::vec & priorTemperature = arma::ones(2));
