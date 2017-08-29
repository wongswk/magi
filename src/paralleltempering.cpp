// [[Rcpp::plugins(cpp11)]]
#include "paralleltempering.h"

std::default_random_engine randgen;
std::uniform_real_distribution<double> unifdistr(0.0,1.0);


void print_info(const arma::umat & swapindicator, const arma::umat & mcmcindicator,
                const vec & temperature, const int & niter) {
  cout << "Parallel tempering finished:\n" 
       << "\tOut of " << niter << " iterations, " 
       << arma::accu(swapindicator.col(0) > 0) << " swap is performed, \n"
       << "swap rate is " << double(arma::accu(swapindicator.col(0) > 0)) / niter
       << endl;
  
  for(int i=0; i<temperature.size(); i++){
    cout << "Chain " << i+1 << " acceptance rate = " 
         << arma::accu(mcmcindicator.row(i))/double(niter) << endl;
  }
  
  cout << "\n========================\n";
  
  for(int i=1; i<temperature.size(); i++){
    int nswap = arma::accu(swapindicator.col(0)==i);
    int nacceptswap = arma::accu(swapindicator.col(0)==i && swapindicator.col(2)==1);
    cout << "Swap between chain " << i << " and chain " << i+1 << ":\n" 
         << "\ttotal swap is " << nswap
         << ", acceptance number = " << nacceptswap
         << ", acceptance rate = " << double(nacceptswap) / double(nswap) << endl;
  }
}

cube parallel_termperingC(std::function<lp (arma::vec)> & lpr, 
                          std::function<mcmcstate (function<lp(vec)>, mcmcstate)> & mcmc, 
                          const arma::vec & temperature, 
                          const arma::vec & initial, 
                          double alpha0, int niter){
  vector<future<mcmcstate>> slave_mcmc(temperature.size());
  vector<future<lp>> slave_eval(temperature.size());
  vector<function<lp(vec)>> lprtempered(temperature.size());
  vector<mcmcstate> paralxs(temperature.size());

  cube retstate(initial.size()+1, temperature.size(), niter);
  arma::umat swapindicator(niter, 3, arma::fill::zeros);
  arma::umat mcmcindicator(temperature.size(), niter, arma::fill::zeros);

  // initial setup
  for(int i=0; i<temperature.size(); i++){
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
  
  for(int i=0; i<temperature.size(); i++){
    paralxs[i].lpv = slave_eval[i].get().value;
  }
  
  for(int it=0; it<niter; it++){
    // MCMC update
    
    for(int i=0; i<temperature.size(); i++){
      slave_mcmc[i] = async(mcmc, lprtempered[i], paralxs[i]);
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
    for(int i=0; i < temperature.size(); i++){
      retstate(0, i, it) = paralxs[i].lpv;
      retstate.slice(it).col(i).subvec(1, initial.size()) = paralxs[i].state;
      mcmcindicator(i, it) = paralxs[i].acc;
    }
  }
  print_info(swapindicator, mcmcindicator, temperature, niter);
  return retstate;
}


mcmcstate metropolis (function<lp(vec)> lpv, mcmcstate current, double stepsize=1.0){
  vec proposal = current.state;
  proposal += arma::randn<vec>(current.state.size())*stepsize;
  
  // cout << proposal << endl;
  
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


// [[Rcpp::export]]
arma::cube main2() {
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
  
  return samples;
}

// [[Rcpp::export]]
arma::cube main3() {
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
  
  return samples;
}

