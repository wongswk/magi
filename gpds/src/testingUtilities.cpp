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
int bandTest(std::string filename="data_band.txt"){
  const int datasize = 201, bandsize = 20;
  
  double xtheta[datasize * 2 + 3];
  double Vmphi[datasize * (bandsize * 2 + 1)];
  double VKinv[datasize * (bandsize * 2 + 1)];
  double VCinv[datasize * (bandsize * 2 + 1)];
  double Rmphi[datasize * (bandsize * 2 + 1)];
  double RKinv[datasize * (bandsize * 2 + 1)];
  double RCinv[datasize * (bandsize * 2 + 1)];
  int bandsizeInput, datasizeInput;
  double cursigma;
  double fnsim[datasize * 2];
  
  ifstream fin(filename);
  for(int i = 0; i < datasize * 2 + 3; i++){
    fin >> xtheta[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> Vmphi[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> VKinv[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> VCinv[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> Rmphi[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> RKinv[i];
  }
  for(int i = 0; i < datasize * (bandsize * 2 + 1); i++){
    fin >> RCinv[i];
  }
  fin >> bandsizeInput;
  fin >> datasizeInput;
  if(bandsizeInput != bandsize || datasizeInput != datasize){
    throw "size not matched";
  }
  fin >> cursigma;
  for(int i = 0; i < datasize * 2; i++){
    fin >> fnsim[i];
  }
  
  double mysum=0;
  cout << std::accumulate(xtheta, xtheta + 405, mysum) << endl;
  
  double ret=0;
  double grad[datasize * 2 + 3];
  
  xthetallik(xtheta, Vmphi, VKinv, VCinv,
             Rmphi, RKinv, RCinv, &bandsizeInput, &datasizeInput,
             &cursigma, fnsim, &ret, grad);
  cout << ret << endl;
  for(int i = 0; i < datasize * 2 + 3; i++){
    cout << grad[i] << ",";
  }
  printf("\ncheck sum = %.20e\n", std::accumulate(grad, grad + datasize * 2 + 3, ret));
  return 0;
};



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

