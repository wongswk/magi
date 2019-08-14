#include <cmath>
#include <random>
#include <vector>
// [[Rcpp::plugins(cpp11)]]
#include <functional>

#include "rcppgpds/classDefinition.h"
#include "rcppgpds/hmc.h"
#include "rcppgpds/tgtdistr.h"
#include "rcppgpds/paralleltempering.h"
#include "rcppgpds/dynamicalSystemModels.h"

gpcov cov_r2cpp(const Rcpp::List & cov_r);
Rcpp::List xthetallikRcpp(const arma::mat & yobs, 
                          const Rcpp::List & covAllDimInput,
                          const arma::vec & sigmaInput,
                          const arma::vec & initial,
                          const std::string modelName,
                          const bool useBand,
                          const Rcpp::NumericVector & priorTemperatureInput);
