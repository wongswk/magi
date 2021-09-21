//
// Created by Shihao Yang on 9/15/21.
//

#define NDEBUG
#define EIGEN_NO_DEBUG

// [[Rcpp::depends(RcppEigen)]]

#include <armadillo>
#include <Eigen/Core>
#include <iostream>
#include <LBFGSB.h>  // Note the different header file

using Eigen::VectorXd;
using namespace LBFGSpp;

class PhiGaussianProcessSmoothing2 {
private:
    int n;
public:
    PhiGaussianProcessSmoothing2(int n_) : n(n_) {}

    double operator()(const Eigen::VectorXd & phisigInput, Eigen::VectorXd & grad){
        grad.fill(1);
        return phisigInput.sum();
    }

};


// [[Rcpp::export]]
Eigen::VectorXd gpsmoothEigen2() {
    const int n = 3;
    // Set up parameters
    LBFGSBParam<double> param;  // New parameter class
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);  // New solver class
    PhiGaussianProcessSmoothing2 fun(n);

    // Bounds
    VectorXd lb = VectorXd::Constant(n, 1e-4);
    VectorXd ub = VectorXd::Constant(n, INFINITY);

    // Initial guess
    VectorXd x = VectorXd::Constant(n, 3.0);

    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return x;
}


// changed function return type and
// the return type of first parameter
double obj_fun_rcpp(double& beta_hat,
                    arma::vec& x, arma::vec& y, arma::vec& tar_val){

    // Changed from % to * as it is only appropriate if
    // `beta_hat` is the same length as x and y.
    // This is because it performs element-wise multiplication
    // instead of a scalar multiplication on a vector
    arma::vec est_val = (pow(x, 2) - pow(y, 2)) * beta_hat;

    // Compute objective value
    double obj_val = sum( pow( est_val - tar_val, 2) );

    // Return a single value
    return obj_val;
}


// [[Rcpp::export]]
arma::vec optim_rcpp(double& init_val,
                     arma::vec& x, arma::vec& y, arma::vec& tar_val){

    // Extract R's optim function
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];

    // Call the optim function from R in C++
    Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
            // Make sure this function is not exported!
                                   Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                   Rcpp::_["method"] = "BFGS",
            // Pass in the other parameters as everything
            // is scoped environmentally
                                   Rcpp::_["x"] = x,
                                   Rcpp::_["y"] = y,
                                   Rcpp::_["tar_val"] = tar_val);

    // Extract out the estimated parameter values
    arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);

    // Return estimated values
    return out;
}
