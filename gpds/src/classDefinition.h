#ifndef CLASSDEFINITION_H
#define CLASSDEFINITION_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <armadillo>

struct hmcstate{
  arma::vec final, finalp, step, trajH;
  double lprvalue, apr, delta;
  int acc;
  arma::mat trajq, trajp;
};

struct lp{
  double value;
  arma::vec gradient;
  lp (const double & tgtv) : value(tgtv) {}
  lp (){}
  lp (const lp & lp2) : value(lp2.value), gradient(lp2.gradient) {}
};

struct mcmcstate {
  arma::vec state;
  double lpv;
  int acc;
  mcmcstate(const hmcstate & x){
    state = x.final;
    lpv = x.lprvalue;
    acc = x.acc;
  }
  mcmcstate(){}
  mcmcstate(const mcmcstate & another){
    state = another.state;
    lpv = another.lpv;
    acc = another.acc;
  }
};

#endif
