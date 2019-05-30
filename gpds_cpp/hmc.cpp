// [[Rcpp::plugins(cpp11)]]
#include "hmc.h"

using namespace arma;

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
//' @noRd
hmcstate basic_hmcC(const std::function<lp (vec)> & lpr, 
                    const vec & initial, 
                    const vec & step, 
                    vec lb, 
                    vec ub,
                    const int nsteps = 1, 
                    const bool traj = false){
  // Check and process the arguments
  if(step.size() != initial.size())
    throw std::runtime_error("step and initial dimention note matched");
  if(nsteps <= 0)
    throw std::runtime_error("Invalid nsteps argument");
  if(lb.size() != initial.size()){
    if(lb.size() == 1){
      double lbval = lb(0);
      lb = vec(initial.size());
      lb.fill(lbval);
    }else{
      throw std::runtime_error("lb and initial dimention note matched");  
    }
  }
  if(ub.size() != initial.size()){
    if(ub.size() == 1){
      double ubval = ub(0);
      ub = vec(initial.size());
      ub.fill(ubval);
    }else{
      throw std::runtime_error("ub and initial dimention note matched");  
    }
  }
  
  
  // Allocate space for the trajectory, if its return is requested.
  mat* trajq = 0;
  mat* trajp = 0;
  vec* trajH = 0;
  if (traj){ 
    trajq = new mat(zeros<mat>(nsteps+1,initial.size()));
    trajp = new mat(zeros<mat>(nsteps+1,initial.size()));
    trajH = new vec(zeros<vec>(nsteps+1));
    // cout << "trajq" << (*trajq)(1,1) << endl;
  }
  
  // Evaluate the log probability and gradient at the initial position
  lp lpx = lpr(initial);
  if(std::isnan(lpx.value)){
    throw std::runtime_error("hmc evaluates the log target density to be NaN at initial value");
  }
  // cout << "Finish Evaluate the log probability and gradient at the initial position" << endl;
  
  // Compute the kinetic energy at the start of the trajectory
  vec initialp = randn<vec>(initial.size());
  double kineticinitial = sum(square(initialp)) / 2.0;
  
  // Compute the trajectory by the leapfrog method
  vec q = initial;
  vec p = initialp;
  vec gr = lpx.gradient;
  if (traj){ 
    (*trajq).row(0) = initial.t();
    (*trajp).row(0) = initialp.t();
    (*trajH)(0) = kineticinitial - lpx.value;
  } 
  double Hinitial = -lpx.value + kineticinitial;
  // cout << "Finish Compute the trajectory by the leapfrog method" << endl;
  
  // Make a half step for momentum at the beginning
  p = p + ((step/2.0) % lpx.gradient);
  
  // Alternate full steps for position and momentum.
  lp lprq;
  for( int i = 0; i < nsteps; i++){
    // Make a full step for the position, and evaluate the gradient at the new position.
    q = q + step % p;
    // Fig 8: Modification to the leapfrog update of q (eq 2.29) to handle constraints
    mat bounces = bouncebyconstraint(q, lb, ub);
    q = bounces.col(0);
    p = p % bounces.col(1);
    
    lprq = lpr(q);
    if(std::isnan(lprq.value)){
      break;
    }
    gr = lprq.gradient;
    
    // Record trajectory if asked to, with half-step for momentum.
    if (traj){ 
      (*trajq).row(i+1) = q.t();
      (*trajp).row(i+1) = (p + (step/2.0) % gr).t();
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
  p = p + (step/2.0) % gr;
  
  // Negate momentum at end of trajectory to make the proposal symmetric.
  p = -p;
  
  // Look at log probability and kinetic energy at the end of the trajectory.
  double lprprop = lprq.value;
  double kineticprop = sum(square(p)) / 2.0;
  
  // Accept or reject the state at the end of the trajectory.
  double Hprop = -lprprop + kineticprop;
  if(std::isnan(Hprop)){
    Hprop = arma::datum::inf;
  }
  
  double delta = Hprop - Hinitial;
  double apr = std::min(1.0,  std::exp(-delta));
  
  // default REJECT 
  vec finalq = initial;
  vec finalp = initialp;
  double lprfinal = lpx.value;
  int acc = 0;
  
  if (as_scalar(randu(1)) < apr) { // ACCEPT
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

lp lpnormal(vec x){
  lp lpx;
  lpx.value = -sum(square(x))/2.0;
  lpx.gradient = -x;
  return lpx;
}





mat bouncebyconstraint(vec x, vec lb, vec ub){
  uvec toolow = x < lb;
  uvec toohigh = x > ub;
  vec uturn = ones<vec>(x.size());
  
  for(unsigned int i = 0; i < toolow.size(); i++){
    if(toolow(i)){
      if(is_finite(ub(i))){
        int k = ceil((lb(i) - x(i))/(2.0*(ub(i)-lb(i))));
        x(i) = x(i) + 2.0*(ub(i)-lb(i))*double(k);
        if(x(i) > ub(i)){
          x(i) = 2.0*ub(i) - x(i);
          uturn(i) = -uturn(i);
        }
      }else{
        x(i) = 2.0*lb(i) - x(i);
        uturn(i) = -uturn(i);
      }
    }
    if(toohigh(i)){
      if(is_finite(lb(i))){
        int k = ceil((lb(i) - x(i))/(2.0*(ub(i)-lb(i))));
        x(i) = x(i) + 2.0*(ub(i)-lb(i))*double(k);
        if(x(i) > ub(i)){
          x(i) = 2.0*ub(i) - x(i);
          uturn(i) = -uturn(i);
        }
      }else{
        x(i) = 2.0*ub(i) - x(i);
        uturn(i) = -uturn(i);
      }
    }
  }
  // cout << x << endl;
  return join_horiz(x,uturn);
}