#include <cmath>
#include <random>
#include <vector>
// [[Rcpp::plugins(cpp11)]]
#include <functional>

#include "classDefinition.h"
#include "hmc.h"
#include "tgtdistr.h"
#include "paralleltempering.h"

gpcov cov_r2cpp(const Rcpp::List & cov_r);