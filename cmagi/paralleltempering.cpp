// [[Rcpp::plugins(cpp11)]]
#include "paralleltempering.h"

using arma::vec;
using arma::mat;
using arma::cube;

void print_info(const arma::umat & swapindicator, const arma::umat & mcmcindicator,
                const vec & temperature, const int & niter) {
  std::cout << "Parallel tempering finished:\n"
       << "\tOut of " << niter << " iterations, " 
       << arma::accu(swapindicator.col(0) > 0) << " swap is performed, \n"
       << "swap rate is " << double(arma::accu(swapindicator.col(0) > 0)) / niter
       << endl;
  
  for(unsigned int i=0; i<temperature.size(); i++){
    std::cout << "Chain " << i+1 << " acceptance rate = "
         << arma::accu(mcmcindicator.row(i))/double(niter) << endl;
  }
  
  std::cout << "\n========================\n";
  
  for(unsigned int i=1; i<temperature.size(); i++){
    int nswap = arma::accu(swapindicator.col(0)==i);
    int nacceptswap = arma::accu(swapindicator.col(0)==i && swapindicator.col(2)==1);
    std::cout << "Swap between chain " << i << " and chain " << i+1 << ":\n"
         << "\ttotal swap is " << nswap
         << ", acceptance number = " << nacceptswap
         << ", acceptance rate = " << double(nacceptswap) / double(nswap) << endl;
  }
}

cube parallel_termperingC(std::function<lp (arma::vec)> & lpr, 
                          std::function<mcmcstate (function<lp(vec)>, mcmcstate)> & mcmc, 
                          const arma::vec & temperature, 
                          const arma::vec & initial, 
                          double alpha0, int niter, bool verbose){
  std::default_random_engine randgen;
  std::uniform_real_distribution<double> unifdistr(0.0,1.0);
  
  vector<future<mcmcstate>> slave_mcmc(temperature.size());
  vector<future<lp>> slave_eval(temperature.size());
  vector<function<lp(vec)>> lprtempered(temperature.size());
  vector<mcmcstate> paralxs(temperature.size());

  cube retstate(initial.size()+1, temperature.size(), niter);
  arma::umat swapindicator(niter, 3, arma::fill::zeros);
  arma::umat mcmcindicator(temperature.size(), niter, arma::fill::zeros);

  // initial setup
  for(unsigned int i=0; i<temperature.size(); i++){
    lprtempered[i] = [&lpr, &temperature, i](vec x) -> lp { 
      lp ret = lpr(x);
      ret.value = ret.value/temperature(i);
      ret.gradient = ret.gradient/temperature(i);
      return ret; 
      };
    paralxs[i].state = initial;
    paralxs[i].acc = 1;
    slave_eval[i] = async(lprtempered[i], paralxs[i].state);
  }
  
  for(unsigned int i=0; i<temperature.size(); i++){
    paralxs[i].lpv = slave_eval[i].get().value;
  }
  
  int nmilestone = niter/10;
  for(int it=0; it<niter; it++){
    if(verbose && it % nmilestone == 0){
      std::cout << "\n === mile stone ===> processed " << it / nmilestone * 10
           << "% of MCMC iterations\n";
    }
    // MCMC update
    
    for(unsigned int i=0; i<temperature.size(); i++){
      // std::cout << paralxs[i].lpv << endl;
      slave_mcmc[i] = async(mcmc, lprtempered[i], paralxs[i]);
    }
    
    for(unsigned int i=0; i<temperature.size(); i++){
      paralxs[i] = slave_mcmc[i].get();
    }
    // swapping
    if(unifdistr(randgen) < alpha0){
      double moveinfo = unifdistr(randgen)*double(temperature.size());
      int movefromid = floor(moveinfo);
      int movetoid = movefromid;
      if(movefromid == 0){
        movetoid ++;
      }else if(movetoid == (int)temperature.size()-1){
        movetoid --;
      }else{
        if(moveinfo - movefromid > 0.5){
          movetoid ++;
        }else{
          movetoid --;
        }
      }
      swapindicator(it, 0) = min(movefromid, movetoid)+1;
      swapindicator(it, 1) = max(movefromid, movetoid)+1;
      
      slave_eval[0] = async(lprtempered[movefromid], paralxs[movetoid].state);
      slave_eval[1] = async(lprtempered[movetoid], paralxs[movefromid].state);
      double log_accp_prob = slave_eval[0].get().value + slave_eval[1].get().value - 
        paralxs[movefromid].lpv - paralxs[movetoid].lpv;
      if(log(unifdistr(randgen)) < log_accp_prob){
        mcmcstate tmp = paralxs[movefromid];
        paralxs[movefromid] = paralxs[movetoid];
        paralxs[movetoid] = tmp;
        swapindicator(it, 2) = 1;
      }
    }
    
    // store states
    for(unsigned int i=0; i < temperature.size(); i++){
      retstate(0, i, it) = paralxs[i].lpv;
      retstate.slice(it).col(i).subvec(1, initial.size()) = paralxs[i].state;
      mcmcindicator(i, it) = paralxs[i].acc;
    }
  }
  if(verbose) {
    print_info(swapindicator, mcmcindicator, temperature, niter); 
  }
  return retstate;
}


mcmcstate metropolis (function<lp(vec)> lpv, mcmcstate current, double stepsize=1.0){
  std::default_random_engine randgen;
  std::uniform_real_distribution<double> unifdistr(0.0,1.0);
  
  vec proposal = current.state;
  proposal += arma::randn<vec>(current.state.size())*stepsize;
  
  // std::cout << proposal << endl;
  
  double proplpv = lpv(proposal).value;
  mcmcstate ret = current;
  ret.acc = 0;
  if(log(unifdistr(randgen)) < proplpv - current.lpv){
    ret.state = proposal;
    ret.lpv = proplpv;
    ret.acc = 1;
  }
  return ret;
}

