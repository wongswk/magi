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
  mat C, Cinv, mphi, Kphi, Kinv, dCdphi1, dCdphi2;
};

gpcov maternCov( vec, mat, int);
gpcov rbfCov( vec, mat, int);
gpcov compact1Cov( vec, mat, int);
lp phisigllik( vec, mat, mat, string kernel = "matern");
lp xthetallik( vec, gpcov, gpcov, double, mat, std::function<mat (vec, mat)> );
lp xthetallik_rescaled( vec, gpcov, gpcov, double, mat, std::function<mat (vec, mat)> );
mat fnmodelODE( vec, mat);