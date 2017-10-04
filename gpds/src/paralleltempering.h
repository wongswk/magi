#include <iostream>
#include <future>
#include <chrono>
#include <random>
#include <armadillo>
#include "hmc.h"
// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>

using namespace std;

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

void print_info(const arma::umat &, const arma::umat &, const arma::vec &, const int &);
arma::cube parallel_termperingC(std::function<lp (arma::vec)> & , 
                          std::function<mcmcstate (function<lp(arma::vec)>, mcmcstate)> &, 
                          const arma::vec &, 
                          const arma::vec &, 
                          double, int);
mcmcstate metropolis (function<lp(arma::vec)>, mcmcstate, double);
  