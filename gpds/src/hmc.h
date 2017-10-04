#ifndef HMC_H
#define HMC_H

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>

using namespace std;
using namespace arma;

struct lp{
  double value;
  vec gradient;
  lp (const double & tgtv) : value(tgtv) {}
  lp (){}
  lp (const lp & lp2) : value(lp2.value), gradient(lp2.gradient) {}
};

struct hmcstate{
  vec final, finalp, step, trajH;
  double lprvalue, apr, delta;
  int acc;
  mat trajq, trajp;
};

hmcstate basic_hmcC(const std::function<lp (vec)> &, 
                    const vec &, 
                    const vec &, 
                    vec, 
                    vec,
                    int, 
                    bool);

lp lpnormal(vec);
mat bouncebyconstraint(vec, vec, vec);
  
#endif