// #include <Rcpp.h>
// using namespace Rcpp;

#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <armadillo>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
int main()
{
  mat A = randu<mat>(4,5);
  mat B = randu<mat>(4,5);
  
  cout << A*B.t() << endl;
  
  return 0;
}

struct lp{
  double value;
  vec gradient;
};

//' basic_hmcC
//' 
//' BASIC HAMILTONIAN MONTE CARLO UPDATE
//' Shihao Yang, 2017, rewriten from Radford M. Neal, 2012.
//'
//' @param lpr       Function returning the log probability of the position part 
//'                  of the state, plus an arbitrary constant, with gradient
//'                  as an attribute if grad=TRUE is passed.
//' @param initial   The initial position part of the state (a vector).
//' @param initialp  The initial momentum part of the state (a vector), default
//'                  is to use momentum variables generated as standard normals.
//' @param nsteps    Number of steps in trajectory used to propose a new state.
//'                  (Default is 1, giving the "Langevin" method.)
//' @param step      Stepsize or stepsizes.  May be scalar or a vector of length 
//'                  equal to the dimensionality of the state.
//' @param traj      TRUE if values of q and p along the trajectory should be 
//'                  returned (default is FALSE).
int basic_hmcC(lp (*lpr)(vec), const vec & initial, vec step,
                int nsteps = 1, bool traj = false){
  // Check and process the arguments
  if(step.size() != initial.size())
    throw "step and initial dimention note matched";
  if(nsteps <= 0)
    throw "Invalid nsteps argument";
  
  // Allocate space for the trajectory, if its return is requested.
  mat* trajq = 0;
  mat* trajp = 0;
  vec* trajH = 0;
  if (traj){ 
    trajq = new mat(zeros<mat>(nsteps+1,initial.size()));
    trajp = new mat(zeros<mat>(nsteps+1,initial.size()));
    trajH = new vec(zeros<vec>(nsteps+1));
    cout << "trajq" << (*trajq)(1,1) << endl;
  }
  
  // Evaluate the log probability and gradient at the initial position
  lp lpx = (*lpr)(initial);
  cout << lpx.value << lpx.gradient << endl;
  
  
  delete trajq;
  delete trajp;
  delete trajH;
  return 0;
}


//' R wrapper for basic_hmcC
// [[Rcpp::export]]
vec hmc(const arma::mat &MgrPos, const arma::cube &OppPos, const arma::mat &ConstitRet){
  return randu<vec>(4);
}

// [[Rcpp::export]]
vec GetMod(vec x, int n){
  vec mod(x.size());
  for(int i=0; i<x.size(); ++i){
    int num = x(i);
    mod(i) = num % n;
  }
  return mod;
}

// [[Rcpp::export]]
vec test(vec x, vec y){
  return x % y;
}

lp lpr(vec x){
  double v = 5.0;
  vec g = zeros<vec>(5);
  lp lpx;
  lpx.value = v;
  lpx.gradient = g;
  return lpx;
}

// [[Rcpp::export]]
int test2(const vec & initial, vec step,
          int nsteps = 1, bool traj = false){
  return basic_hmcC(lpr, initial, step, nsteps, traj);
}
