#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
// #include <armadillo>
// [[Rcpp::depends(RcppArmadillo)]]
#include <armadillo>

#include "classDefinition.h"

using namespace std;

gpcov maternCov( arma::vec, arma::mat, int);
gpcov rbfCov( arma::vec, arma::mat, int);
gpcov compact1Cov( arma::vec, arma::mat, int);
lp phisigllik( arma::vec, arma::mat, arma::mat, string kernel = "matern");
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
arma::mat fnmodelODE(const arma::vec &, const arma::mat &);