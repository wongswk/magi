// [[Rcpp::plugins(cpp11)]]
// #include "paralleltempering.h"
#include <iostream>
#include <future>
#include <chrono>
#include <random>
#include <armadillo>
// #include "hmc.h"

using namespace std;
using arma::vec;
using arma::mat;
using arma::cube;

std::default_random_engine randgen;
std::uniform_real_distribution<double> unifdistr(0.0,1.0);

int square(int x) {
  return x * x;
}


struct mcmcstate {
  vec state;
  double lpv;
};

mat parallel_termperingC(std::function<double (arma::vec)> & lpv, 
                          std::function<mcmcstate (function<double(vec)>, mcmcstate)> & mcmc, 
                          const arma::vec & temperature, 
                          const arma::vec & initial, 
                          double alpha0, int niter){
  vector<future<mcmcstate>> slave_mcmc(temperature.size());
  vector<future<double>> slave_eval(temperature.size());
  vector<function<double(vec)>> lpvtempered(temperature.size());
  vector<mcmcstate> paralxs(temperature.size());

  cube retstate(initial.size()+1, temperature.size(), niter);
  
  cout << "begin initial setup" << endl;
  // initial setup
  for(int i=0; i<temperature.size(); i++){
    cout << "i = " << i << endl;
    lpvtempered[i] = [&lpv, &temperature, i](vec x) -> double { return lpv(x)/temperature(i); };
    paralxs[i].state = initial;
    slave_eval[i] = async(lpvtempered[i], paralxs[i].state);
  }
  
  cout << "finish slave" << endl;
  cout << "length of slave_eval = " << slave_eval.size() << endl;
  
  for(int i=0; i<temperature.size(); i++){
    paralxs[i].lpv = slave_eval[i].get();
    cout << "paralxs["<< i << "].lpv = " << paralxs[i].lpv << endl;
  }
  
  cout << "finish initial setup" << endl;
  
  
  
  for(int it=0; it<niter; it++){
    // MCMC update
    
    cout << "paralxs[0] =\n" << paralxs[0].lpv << paralxs[0].state << endl;
    cout << "mcmc eval =\n" << mcmc(lpvtempered[0], paralxs[0]).state << endl;
    cout << "finish mcmc eval" << endl;
    
    for(int i=0; i<temperature.size(); i++){
      cout << i << endl;
      slave_mcmc[i] = async(mcmc, lpvtempered[i], paralxs[i]);
    }
    cout << "finish async setup" << endl;
    
    for(int i=0; i<temperature.size(); i++){
      paralxs[i] = slave_mcmc[i].get();
    }
    cout << "finish MCMC update" << endl;
    // swapping
    if(unifdistr(randgen) < alpha0){
      double moveinfo = unifdistr(randgen)*double(temperature.size());
      int movefromid = floor(moveinfo);
      int movetoid = movefromid;
      if(movefromid == 0){
        movetoid ++;
      }else if(movetoid == temperature.size()-1){
        movetoid --;
      }else{
        if(moveinfo - movefromid > 0.5){
          movetoid ++;
        }else{
          movetoid --;
        }
      }
      slave_eval[0] = async(lpvtempered[movefromid], paralxs[movetoid].state);
      slave_eval[1] = async(lpvtempered[movetoid], paralxs[movefromid].state);
      double log_accp_prob = slave_eval[0].get() + slave_eval[1].get() - 
        paralxs[movefromid].lpv - paralxs[movetoid].lpv;
      if(log(unifdistr(randgen)) < log_accp_prob){
        mcmcstate tmp = paralxs[movefromid];
        paralxs[movefromid] = paralxs[movetoid];
        paralxs[movetoid] = tmp;
      }
    }
    
    // store states
    for(int i=0; i < temperature.size(); i++){
      retstate(0, i, it) = paralxs[i].lpv;
      retstate.slice(it).col(i).subvec(1, initial.size()) = paralxs[i].state;
    }
  }
  return retstate;
}


mcmcstate metropolis (function<double(vec)> lpv, mcmcstate current, double stepsize=1.0){
  vec proposal = current.state;
  proposal += arma::randn<vec>(current.state.size())*stepsize;
  
  cout << proposal << endl;
  
  double proplpv = lpv(proposal);
  mcmcstate ret = current;
  if(log(unifdistr(randgen)) < proplpv - current.lpv){
    ret.state = proposal;
    ret.lpv = proplpv;
  }
  return ret;
}

// [[Rcpp::export]]
arma::mat main2() {
  function<double(vec)> lpnormalvalue = [](vec x) {return -arma::sum(arma::square(x))/2.0;};
  vec temperature = arma::linspace<vec>(8, 1, 8);
  function<mcmcstate(function<double(vec)>, mcmcstate)> metropolis_tuned =
    std::bind(metropolis, std::placeholders::_1, std::placeholders::_2, 1.0);
  
  mat samples = parallel_termperingC(lpnormalvalue, 
                                     metropolis_tuned, 
                                     temperature, 
                                     arma::zeros<vec>(4), 
                                     0.05, 
                                     1e4);
  
  return samples;
}

arma::mat main3() {
  function<double(vec)> lpnormalvalue = [](vec x) {return -arma::sum(arma::square(x))/2.0;};
  vec temperature = arma::linspace<vec>(8, 1, 8);
  function<mcmcstate(function<double(vec)>, mcmcstate)> metropolis_tuned =
    std::bind(metropolis, std::placeholders::_1, std::placeholders::_2, 1.0);
  
  mcmcstate initial;
  initial.state = arma::zeros<vec>(4);
  initial.lpv = lpnormalvalue(initial.state);
  return metropolis_tuned(lpnormalvalue, initial).state;
}

int main() {
  cout << main2() << endl;
  return 0;
}
