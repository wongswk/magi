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

struct hmcstate{
  vec final, finalp, step, trajH;
  double lprvalue, apr, delta;
  int acc;
  mat trajq, trajp;
  lp lprfinal;
};

const std::default_random_engine randgen;
std::uniform_real_distribution<double> unifdistr(0.0,1.0);

//' basic_hmcC
//' 
//' BASIC HAMILTONIAN MONTE CARLO UPDATE
//' Shihao Yang, 2017, rewriten from Radford M. Neal, 2012.
//'
//' @param lpr       Function returning the log probability of the position part 
//'                  of the state, plus an arbitrary constant, with gradient
//'                  as an attribute if grad=TRUE is passed.
//' @param initial   The initial position part of the state (a vector).
//' @param nsteps    Number of steps in trajectory used to propose a new state.
//'                  (Default is 1, giving the "Langevin" method.)
//' @param step      Stepsize or stepsizes.  May be scalar or a vector of length 
//'                  equal to the dimensionality of the state.
//' @param traj      TRUE if values of q and p along the trajectory should be 
//'                  returned (default is FALSE).
hmcstate basic_hmcC(lp (*lpr)(vec), const vec & initial, vec step,
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
  
  // Compute the kinetic energy at the start of the trajectory
  vec initialp = randn<vec>(initial.size());
  double kineticinitial = sum(square(initialp)) / 2.0;
  
  // Compute the trajectory by the leapfrog method
  vec q = initial;
  vec p = initialp;
  vec gr = lpx.gradient;
  if (traj){ 
    (*trajq).row(0) = initial;
    (*trajp).row(0) = initialp;
    (*trajH)(0) = kineticinitial - lpx.value;
  } 
  double Hinitial = -lpx.value + kineticinitial;
  
  // Make a half step for momentum at the beginning
  p = p + ((step/2.0) % lpx.gradient);
  
  // Alternate full steps for position and momentum.
  lp lprq;
  for( int i = 0; i < nsteps; i++){
    // Make a full step for the position, and evaluate the gradient at the new position.
    q = q + step % p;
    lprq = lpr(q);
    gr = lprq.gradient;
    
    // Record trajectory if asked to, with half-step for momentum.
    if (traj){ 
      (*trajq).row(i+1) = q;
      (*trajp).row(i+1) = p + (step/2.0) % gr;
      (*trajH)(i+1) = sum(square( (*trajp).row(i+1) ))/2.0 - lprq.value;
      
      if ((*trajH)(i+1) - Hinitial > 50.0) {
        break;
      }
    }
    
    // Make a full step for the momentum, except when we're coming to the end of the trajectory.  
    if (i != nsteps-1)
      p = p + step % gr;
  }
  // Make a half step for momentum at the end.  
  p = p + (step/2.0) * gr;
  
  // Negate momentum at end of trajectory to make the proposal symmetric.
  p = -p;
  
  // Look at log probability and kinetic energy at the end of the trajectory.
  double lprprop = lprq.value;
  double kineticprop = sum(square(p)) / 2.0;
  
  // Accept or reject the state at the end of the trajectory.
  double Hprop = -lprprop + kineticprop;
  double delta = Hprop - Hinitial;
  double apr = std::min(1.0,  std::exp(-delta));
  
  // default REJECT 
  vec finalq = initial;
  vec finalp = initialp;
  double lprfinal = lpx.value;
  int acc = 0;
  
  if (double(unifdistr(randgen)) < apr) { // ACCEPT
    finalq = q;
    finalp = p;
    lprfinal = lprq.value;
    acc = 1;
  }
  
  // Return new state, its log probability and gradient, plus additional
  // information, including the trajectory, if requested.
  
  hmcstate ret;
  ret.final = finalq;
  ret.finalp = finalp;
  ret.lprvalue = lprfinal;
  ret.step = step;
  ret.apr = apr;
  ret.acc = acc;
  ret.delta = delta;
  
  if (traj) { 
    ret.trajq = (*trajq);
    ret.trajp = (*trajp);
    ret.trajH = (*trajH);
  }
  
  // clean up and return;
  delete trajq;
  delete trajp;
  delete trajH;
  
  return ret;
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

lp lpnormal(vec x){
  lp lpx;
  lpx.value = sum(square(x))/2.0;
  lpx.gradient = x;
  return lpx;
}

//' R wrapper for basic_hmcC
// [[Rcpp::export]]
Rcpp::List hmc(const vec & initial, vec step,
          int nsteps = 1, bool traj = false){
  hmcstate post = basic_hmcC(lpnormal, initial, step, nsteps, traj);
  
  return List::create(Named("final")=post.final,
                      Named("final.p")=post.finalp,
                      Named("lpr")=post.lprvalue,
                      Named("step")=post.step,
                      Named("apr")=post.apr,
                      Named("acc")=post.acc,
                      Named("delta")=post.delta);
}
