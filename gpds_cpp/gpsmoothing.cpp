#include <armadillo>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <cppoptlib/solver/bfgssolver.h>

#include "tgtdistr.h"


// [[Rcpp::export]]
arma::vec calcFrequencyBasedPrior(const arma::vec & x){
    arma::cx_mat z = arma::fft(x);
    arma::vec zmod(z.n_rows);
    for(unsigned int i = 0; i < z.n_rows; i++){
        zmod(i) = std::norm(z(i, 0));
    }
    const arma::vec & zmodEffective = zmod(arma::span(1, (zmod.size() - 1) / 2));
    const arma::vec zmodEffectiveSorted = arma::sort(zmodEffective);
    double upperQuarter = zmodEffectiveSorted(static_cast<unsigned int>(std::ceil(zmodEffectiveSorted.size() * 0.75)));
    double lowerQuarter = zmodEffectiveSorted(static_cast<unsigned int>(std::floor(zmodEffectiveSorted.size() * 0.25)));
    double iqr = upperQuarter - lowerQuarter;
    arma::vec outliers = zmodEffective(zmodEffective > upperQuarter + 1.5 * iqr);
    long long int freq = 1;
    if(!outliers.empty()){
        for(unsigned long long int i = zmodEffective.size(); i > 0; i--){
            if(arma::min(arma::abs(zmodEffective(i - 1) - outliers)) < 1e-6){
                freq = i;
                break;
            }
        }

    }else{
        freq = zmodEffective.index_max();
    }

    double meanFactor = 0.5 / freq;
    double sdFactor = (1 - meanFactor) / 3;

    return arma::vec({meanFactor, sdFactor});
}


class PhiGaussianProcessSmoothing : public cppoptlib::BoundedProblem<double> {
public:
    std::string kernel;
    const arma::mat & yobs;
    const arma::mat & dist;
    const unsigned int numparam;
    const bool useFrequencyBasedPrior;
    arma::vec priorFactor;
    double maxDist;

    double value(const Eigen::VectorXd & phisigInput) override {
        if ((phisigInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        const arma::vec & phisig = arma::vec(const_cast<double*>(phisigInput.data()), numparam, false, false);
        const lp & out = phisigllik(phisig, yobs, dist, kernel);
        double penalty = 0;
        if (useFrequencyBasedPrior) {
            for (unsigned j = 0; j < yobs.n_cols; j++){
                penalty += -0.5 * std::pow((phisig(2*j+1) - maxDist * priorFactor(0)) / (maxDist * priorFactor(1)), 2);
            }
        }
        return -(out.value + penalty);
    }

    void gradient(const Eigen::VectorXd & phisigInput, Eigen::VectorXd & grad) override {
        if ((phisigInput.array() < this->lowerBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < numparam; i++){
                if(phisigInput[i] < this->lowerBound()[i]){
                    grad[i] = -1;
                }
            }
            return;
        }
        const arma::vec & phisig = arma::vec(const_cast<double*>(phisigInput.data()), numparam, false, false);
        const lp & out = phisigllik(phisig, yobs, dist, kernel);
        for(unsigned i = 0; i < numparam; i++){
            grad[i] = -out.gradient(i);
        }
        double penalty = 0;
        if (useFrequencyBasedPrior) {
            for (unsigned j = 0; j < yobs.n_cols; j++){
                penalty = (phisig(2*j+1) - maxDist * priorFactor(0)) / std::pow((maxDist * priorFactor(1)), 2);
                grad[2*j+1] += penalty;
            }
        }
    }

    PhiGaussianProcessSmoothing(const arma::mat & yobsInput,
                                const arma::mat & distInput,
                                std::string kernelInput,
                                const unsigned int numparamInput,
                                const bool useFrequencyBasedPriorInput) :
            BoundedProblem(numparamInput),
            kernel(std::move(kernelInput)),
            yobs(yobsInput),
            dist(distInput),
            numparam(numparamInput),
            useFrequencyBasedPrior(useFrequencyBasedPriorInput) {
        Eigen::VectorXd lb(numparam);
        lb.fill(1e-4);
        this->setLowerBound(lb);
        maxDist = dist.max();

        priorFactor = arma::zeros(2);
        if(useFrequencyBasedPrior){
            for (unsigned j = 0; j < yobs.n_cols; j++){
                priorFactor += calcFrequencyBasedPrior(yobs);
            }
            priorFactor /= yobs.n_cols;
        }
    }
};


// [[Rcpp::export]]
arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   const unsigned int numparamInput,
                   bool useFrequencyBasedPrior = false) {
    PhiGaussianProcessSmoothing objective(yobsInput, distInput, std::move(kernelInput), numparamInput, useFrequencyBasedPrior);
    cppoptlib::LbfgsbSolver<PhiGaussianProcessSmoothing> solver;
    Eigen::VectorXd phisig(2 * yobsInput.n_cols + 1);
    double maxDist = distInput.max();
    double sdOverall = 0;
    for(unsigned i = 0; i < yobsInput.n_cols; i++) {
        phisig[2 * i] = 0.5 * arma::stddev(yobsInput.col(i));
        phisig[2 * i + 1] = 0.5 * maxDist;
        sdOverall += phisig[2 * i];
    }
    phisig[2 * yobsInput.n_cols] = sdOverall / yobsInput.n_cols;
    solver.minimize(objective, phisig);
    const arma::vec & phisigArgmin = arma::vec(phisig.data(), 2 * yobsInput.n_cols + 1, false, false);
    return phisigArgmin;
}
