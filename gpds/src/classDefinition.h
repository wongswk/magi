#ifndef CLASSDEFINITION_H
#define CLASSDEFINITION_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// #include <armadillo> // uncomment this for pure c++ compilation 
#include <iostream>
#include <stdio.h>

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

struct gpcov {
  arma::mat C, Cinv, mphi, Kphi, Kinv, CeigenVec, KeigenVec, mphiLeftHalf;
  arma::cube dCdphiCube;
  arma::mat CinvBand, mphiBand, KinvBand;
  arma::vec Ceigen1over, Keigen1over, mu, dotmu;
  int bandsize;
};

class OdeSystem {
public:
  // row is observations, col is each X variable
  std::function<arma::mat (arma::vec, arma::mat)> fOde;
  // row is observations, col is each partial X denominator, slice is each X variable numerator
  std::function<arma::cube (arma::vec, arma::mat)> fOdeDx;
  // row is observations, col is each partial theta denominator, slice is each X variable numerator
  std::function<arma::cube (arma::vec, arma::mat)> fOdeDtheta;
  std::string name;
  
  OdeSystem(
    const std::function<arma::mat (arma::vec, arma::mat)> & fOdeInput,
    const std::function<arma::cube (arma::vec, arma::mat)> & fOdeDxInput,
    const std::function<arma::cube (arma::vec, arma::mat)> & fOdeDthetaInput
  ) : fOde(fOdeInput), fOdeDx(fOdeDxInput), fOdeDtheta(fOdeDthetaInput) {};
  OdeSystem() {};
};

#endif
