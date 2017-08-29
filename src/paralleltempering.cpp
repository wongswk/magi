// [[Rcpp::plugins(cpp11)]]
// #include "paralleltempering.h"
#include <iostream>
#include <future>
#include <chrono>
#include <random>
#include <armadillo>

using namespace std;
using arma::vec;
using arma::mat;
using arma::cube;

std::default_random_engine randgen;
std::uniform_real_distribution<double> unifdistr(0.0,1.0);

int square(int x) {
  return x * x;
}

int main() {
  vector<future<int>> a;
  for(int i=0; i<20; i++){
    a.push_back(async(&square, i));
  }
  int v=0;
  cout << "Main thread " << endl;
  for(auto& aeach : a){
    v += aeach.get();
  }
  cout << "The thread returned " << v << endl;
  return 0;
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
  vector<future<mcmcstate>> slave_mcmc;
  vector<future<double>> slave_eval;
  vector<function<double(vec)>> lpvtempered;
  vector<mcmcstate> paralxs;

  cube retstate(initial.size()+1, temperature.size(), niter);
  
  // initial setup
  for(int i=0; i<temperature.size(); i++){
    lpvtempered[i] = [&](vec x) -> double { return lpv(x)/temperature(i); };
    paralxs[i].state = initial;
    slave_eval.push_back(async(lpvtempered[i], paralxs[i].state));
  }
  
  for(int i=0; i<temperature.size(); i++){
    paralxs[i].lpv = slave_eval[i].get();
  }
  
  for(int it=0; it<niter; it++){
    // MCMC update
    for(int i=0; i<temperature.size(); i++){
      slave_mcmc.push_back(async(mcmc, lpvtempered[i], paralxs[i]));
    }
    
    for(int i=0; i<temperature.size(); i++){
      paralxs[i] = slave_mcmc[i].get();
    }
    
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
      if(unifdistr(randgen) < log_accp_prob){
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