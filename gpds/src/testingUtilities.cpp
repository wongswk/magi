#include "hmc.h"
#include "band.h"
#include "paralleltempering.h"

using namespace arma;

// [[Rcpp::export]]
int hmcTest(){
  arma::vec initial = arma::zeros<arma::vec>(4);
  arma::vec step(4);
  step.fill(0.05);
  int nsteps = 20;
  bool traj = true;
  hmcstate post = basic_hmcC(lpnormal, initial, step, 
                             {-arma::datum::inf}, 
                             {arma::datum::inf}, 
                             nsteps, traj);
  // for(int i; i < post.final.size(); i++)
    //   cout << post.final(i) << endl;
  cout << post.final << endl;
  return 0;
}

// [[Rcpp::export]]
int bandTest(){
  std::cout << "in bandTest\n";
  return mainBand();
}



// [[Rcpp::export]]
arma::cube paralleltemperingTest1() {
  function<lp(vec)> lpnormal = [](vec x) {return lp(-arma::sum(arma::square(x))/2.0);};
  vec temperature = arma::linspace<vec>(8, 1, 8);
  std::function<mcmcstate(function<lp(vec)>, mcmcstate)> metropolis_tuned =
    std::bind(metropolis, std::placeholders::_1, std::placeholders::_2, 1.0);
  
  cube samples = parallel_termperingC(lpnormal, 
                                      metropolis_tuned, 
                                      temperature, 
                                      arma::zeros<vec>(4), 
                                      0.05, 
                                      1e4);
  // FIXME segfault 'memory not mapped' when use with the package, 
  // but no error if simply sourcecpp
  return samples;
}

// [[Rcpp::export]]
arma::cube paralleltemperingTest2() {
  function<lp(vec)> lpnormalvalue = [](vec x) {
    return lp(log(exp(-arma::sum(arma::square(x+4))/2.0) + exp(-arma::sum(arma::square(x-4))/2.0)));
  };
  vec temperature = {1, 1.3, 1.8, 2.5, 3.8, 5.7, 8};
  function<mcmcstate(function<lp(vec)>, mcmcstate)> metropolis_tuned =
    std::bind(metropolis, std::placeholders::_1, std::placeholders::_2, 1.0);
  
  cube samples = parallel_termperingC(lpnormalvalue, 
                                      metropolis_tuned, 
                                      temperature, 
                                      arma::zeros<vec>(4), 
                                      0.125, 
                                      1e5);
  cout << "parallel_termperingC finished, before returning from paralleltemperingTest2\n";
  // FIXME segfault 'memory not mapped' when use with the package, 
  // but no error if simply sourcecpp
  return samples;
}

