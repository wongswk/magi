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
               const gpcov & CovV,
               const gpcov & CovR,
               const double & sigma,
               const arma::mat & yobs,
               const std::function<arma::mat (arma::vec, arma::mat)> & fODE);
lp xthetallik_rescaled( const arma::vec & xtheta,
                        const gpcov & CovV,
                        const gpcov & CovR,
                        const double & sigma,
                        const arma::mat & yobs,
                        const std::function<arma::mat (arma::vec, arma::mat)> & fODE);
lp xthetallikBandApprox( arma::vec xtheta, 
                         gpcov & CovV, 
                         gpcov & CovR, 
                         double & sigma, 
                         arma::mat & yobs);
arma::mat fnmodelODE(const arma::vec &, const arma::mat &);
