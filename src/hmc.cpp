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
// [[Rcpp::export]]
void basic_hmcC(double (*lpr)(vec), const vec & initial, int nsteps = 1,
                double step = 1.0, bool traj = false){
  
}