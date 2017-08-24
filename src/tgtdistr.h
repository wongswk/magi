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

gpcov maternCov( vec, mat, int);
lp phisigllik( vec, mat, mat);