#ifndef HMC_H
#define HMC_H

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>
#include "classDefinition.h"

using namespace std;

hmcstate basic_hmcC(const std::function<lp (arma::vec)> &,
                    const arma::vec &,
                    const arma::vec &,
                    arma::vec,
                    arma::vec,
                    int,
                    bool);

lp lpnormal(arma::vec);
arma::mat bouncebyconstraint(arma::vec, arma::vec, arma::vec);

#endif
