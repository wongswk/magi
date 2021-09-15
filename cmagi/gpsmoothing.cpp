#define NDEBUG
#define EIGEN_NO_DEBUG

#include <armadillo>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <Eigen/Core>

// [[Rcpp::depends(RcppEigen)]]

class PhiGaussianProcessSmoothing : public cppoptlib::BoundedProblem<double> {
public:
    const unsigned int numparam;

    double value(const Eigen::VectorXd & phisigInput) override {
        return phisigInput.sum();
    }

    void gradient(const Eigen::VectorXd & phisigInput, Eigen::VectorXd & grad) override {
        grad.fill(1);
        return;
    }

    PhiGaussianProcessSmoothing(const unsigned int numparamInput) :
            BoundedProblem(numparamInput),
            numparam(numparamInput) {
        Eigen::VectorXd lb(numparam);
        lb.fill(1e-4);
        this->setLowerBound(lb);

        Eigen::VectorXd ub(numparam);
        ub.fill(10);
    }
};


// [[Rcpp::export]]
arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   const double sigmaExogenScalar = -1,
                   bool useFrequencyBasedPrior = false) {
    unsigned int phiDim;
    if(kernelInput == "generalMatern") {
        phiDim = 2;
    }else if(kernelInput == "matern") {
        phiDim = 2;
    }else if(kernelInput == "compact1") {
        phiDim = 2;
    }else if(kernelInput == "periodicMatern"){
        phiDim = 3;
    }else{
        throw std::invalid_argument("kernelInput invalid");
    }
    unsigned int numparam;

    if(sigmaExogenScalar > 0){
        numparam = phiDim * yobsInput.n_cols;
    }else{
        numparam = phiDim * yobsInput.n_cols + 1;
    }

    PhiGaussianProcessSmoothing objective(numparam);
    cppoptlib::LbfgsbSolver<PhiGaussianProcessSmoothing> solver;
    // phi sigma 1st initial value for optimization
    Eigen::VectorXd phisigAttempt1(numparam);
    phisigAttempt1.fill(1);
    solver.minimize(objective, phisigAttempt1);

    arma::vec phisigArgmin(numparam);
    for(unsigned i = 0; i < numparam; i++){
        phisigArgmin(i) = phisigAttempt1[i];
    }
    return phisigArgmin;
}

// [[Rcpp::export]]
Eigen::VectorXd gpsmoothEigen() {
    unsigned int numparam = 3;

    PhiGaussianProcessSmoothing objective(numparam);
    cppoptlib::LbfgsbSolver<PhiGaussianProcessSmoothing> solver;
    // phi sigma 1st initial value for optimization
    Eigen::VectorXd phisigAttempt1(numparam);
    phisigAttempt1.fill(1);
    solver.minimize(objective, phisigAttempt1);

    return phisigAttempt1;
}
