//
// Created by Shihao Yang on 9/15/21.
//

#define NDEBUG
#define EIGEN_NO_DEBUG

// [[Rcpp::depends(RcppEigen)]]

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
