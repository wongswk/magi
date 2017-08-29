#include <iostream>
#include <future>
#include <chrono>
#include <random>
#include <armadillo>
#include "hmc.h"
// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>

using namespace std;
using arma::vec;
using arma::mat;
using arma::cube;

struct mcmcstate {
  vec state;
  double lpv;
  int acc;
  mcmcstate(const hmcstate & x){
    state = x.final;
    lpv = x.lprvalue;
    acc = x.acc;
  }
};

void print_info(const arma::umat &, const arma::umat &, const vec &, const int &);
cube parallel_termperingC(std::function<double (arma::vec)> & lpv, 
                          std::function<mcmcstate (function<double(vec)>, mcmcstate)> & mcmc, 
                          const arma::vec & temperature, 
                          const arma::vec & initial, 
                          double alpha0, int niter);
mcmcstate metropolis (function<double(vec)>, mcmcstate, double);
arma::cube main2();
arma::cube main3();
  