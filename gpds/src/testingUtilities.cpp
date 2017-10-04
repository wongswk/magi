#include "hmc.h"

// [[Rcpp::export]]
int hmcTest(){
  vec initial = zeros<vec>(4);
  vec step(4);
  step.fill(0.05);
  int nsteps = 20;
  bool traj = true;
  hmcstate post = basic_hmcC(lpnormal, initial, step, 
                             {-datum::inf}, 
                             {datum::inf}, 
                             nsteps, traj);
  // for(int i; i < post.final.size(); i++)
    //   cout << post.final(i) << endl;
  cout << post.final << endl;
  return 0;
}


