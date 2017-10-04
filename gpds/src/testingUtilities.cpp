#include "hmc.h"
#include "band.h"

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
