#ifndef HMC_H
#define HMC_H

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>

using namespace std;

struct lp{
  double value;
  arma::vec gradient;
  lp (const double & tgtv) : value(tgtv) {}
  lp (){}
  lp (const lp & lp2) : value(lp2.value), gradient(lp2.gradient) {}
};

struct hmcstate{
  arma::vec final, finalp, step, trajH;
  double lprvalue, apr, delta;
  int acc;
  arma::mat trajq, trajp;
};

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
