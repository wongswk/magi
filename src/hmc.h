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
};

struct hmcstate{
  vec final, finalp, step, trajH;
  double lprvalue, apr, delta;
  int acc;
  mat trajq, trajp;
  lp lprfinal;
};

hmcstate basic_hmcC(lp (*lpr)(vec), const vec &, vec, int, bool);
int main();
lp lpnormal(vec);

#endif