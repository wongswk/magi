//
// Created by Shihao Yang on 5/30/19.
//
#include "../hmc.h"
#include "../band.h"
#include "../paralleltempering.h"
#include "../tgtdistr.h"
#include "../dynamicalSystemModels.h"

using namespace arma;

int main(){
    arma::vec initial = arma::zeros<arma::vec>(4);
    arma::vec step(4);
    step.fill(0.05);
    int nsteps = 20;
    bool traj = true;
//    hmcstate post = basic_hmcC(lpnormal, initial, step,
//                               {-arma::datum::inf},
//                               {arma::datum::inf},
//                               nsteps, traj);
    // for(int i; i < post.final.size(); i++)`
    //   std::cout << post.final(i) << endl;
    // std::cout << post.final << endl;
    return 0;
}
