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

class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};


// [[Rcpp::export]]
Eigen::VectorXd gpsmoothEigen2() {
    const int n = 10;
    // Set up parameters
    LBFGSBParam<double> param;  // New parameter class
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);  // New solver class
    Rosenbrock fun(n);

    // Bounds
    VectorXd lb = VectorXd::Constant(n, 2.0);
    VectorXd ub = VectorXd::Constant(n, 4.0);

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
