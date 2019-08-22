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


// [[Rcpp::export]]
arma::cube calcMeanCurve(const arma::vec & xInput,
                        const arma::vec & yInput,
                        const arma::vec & xOutput,
                        const arma::mat & phiCandidates,
                        const arma::vec & sigmaCandidates,
                        const std::string kerneltype = "generalMatern",
                        const bool useDeriv = false) {
    assert(kerneltype == "generalMatern");

    const arma::vec & tvec = arma::join_vert(xOutput, xInput);
    arma::mat distSigned(tvec.size(), tvec.size());
    for(unsigned int i = 0; i < distSigned.n_cols; i++){
        distSigned.col(i) = tvec - tvec(i);
    }
    arma::cube ydyOutput(xOutput.size(), phiCandidates.n_cols, 2, arma::fill::zeros);
    arma::mat & yOutput = ydyOutput.slice(0);
    arma::mat & dyOutput = ydyOutput.slice(1);

    int complexity = 0;
    if(useDeriv){
        complexity = 3;
    }

    for(unsigned it = 0; it < phiCandidates.n_cols; it++){
        const double & sigma = sigmaCandidates(it);
        const arma::vec & phi = phiCandidates.col(it);

        gpcov covObj = generalMaternCov(phi, distSigned, complexity);
        arma::mat C = std::move(covObj.C);

        arma::vec Cdiag = C.diag();
        Cdiag(arma::span(xOutput.size(), Cdiag.size() - 1)) += std::pow(sigma, 2);
        C.diag() = Cdiag;
        yOutput.col(it) = C(arma::span(0, xOutput.size() - 1),
                            arma::span(xOutput.size(), Cdiag.size() - 1)) *
                          arma::solve(C(arma::span(xOutput.size(), Cdiag.size() - 1),
                                        arma::span(xOutput.size(), Cdiag.size() - 1)),
                                      yInput);
        if(useDeriv){
            dyOutput.col(it) = covObj.Cprime(arma::span(0, xOutput.size() - 1),
                                             arma::span(0, xOutput.size() - 1)) *
                               arma::solve(C(arma::span(0, xOutput.size() - 1),
                                             arma::span(0, xOutput.size() - 1)),
                                           yOutput.col(it));
        }
    }
    return ydyOutput;
}


class ThetaOptim : public cppoptlib::BoundedProblem<double> {
public:
    const arma::mat & yobs;
    const OdeSystem & fOdeModel;
    const std::vector<gpcov> & covAllDimensions;
    const arma::vec & sigmaAllDimensions;
    const arma::vec & priorTemperature;
    const arma::mat & xInit;
    const unsigned int numparam;

    double value(const Eigen::VectorXd & thetaInput) override {
        if ((thetaInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        if ((thetaInput.array() > this->upperBound().array()).any()){
            return INFINITY;
        }
        const arma::vec & xtheta = arma::join_vert(
            arma::vectorise(xInit),
            arma::vec(const_cast<double*>(thetaInput.data()), numparam, false, false)
        );
        const lp & out = xthetallik(
                xtheta,
                covAllDimensions,
                sigmaAllDimensions,
                yobs,
                fOdeModel,
                true,
                priorTemperature);
        return -out.value;
    }

    void gradient(const Eigen::VectorXd & thetaInput, Eigen::VectorXd & grad) override {
        if ((thetaInput.array() < this->lowerBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < numparam; i++){
                if(thetaInput[i] < this->lowerBound()[i]){
                    grad[i] = -1;
                }
            }
            return;
        }
        if ((thetaInput.array() > this->upperBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < numparam; i++){
                if(thetaInput[i] > this->upperBound()[i]){
                    grad[i] = 1;
                }
            }
            return;
        }
        const arma::vec & xtheta = arma::join_vert(
                arma::vectorise(xInit),
                arma::vec(const_cast<double*>(thetaInput.data()), numparam, false, false)
        );
        const lp & out = xthetallik(
                xtheta,
                covAllDimensions,
                sigmaAllDimensions,
                yobs,
                fOdeModel,
                true,
                priorTemperature);
        for(unsigned i = 0; i < numparam; i++){
            grad[i] = -out.gradient(i);
        }
    }

    ThetaOptim(const arma::mat & yobsInput,
               const OdeSystem & fOdeModelInput,
               const std::vector<gpcov> & covAllDimensionsInput,
               const arma::vec & sigmaAllDimensionsInput,
               const arma::vec & priorTemperatureInput,
               const arma::mat & xInitInput,
               const unsigned int numparamInput) :
            BoundedProblem(numparamInput),
            yobs(yobsInput),
            fOdeModel(fOdeModelInput),
            covAllDimensions(covAllDimensionsInput),
            sigmaAllDimensions(sigmaAllDimensionsInput),
            priorTemperature(priorTemperatureInput),
            xInit(xInitInput),
            numparam(numparamInput) {
        const Eigen::Map<Eigen::VectorXd> lb (const_cast<double*>(fOdeModel.thetaLowerBound.memptr()), numparam);
        this->setLowerBound(lb);
        const Eigen::Map<Eigen::VectorXd> ub (const_cast<double*>(fOdeModel.thetaUpperBound.memptr()), numparam);
        this->setUpperBound(ub);
    }
};

// [[Rcpp::export]]
arma::vec optimizeThetaInit(const arma::mat & yobsInput,
                            const OdeSystem & fOdeModelInput,
                            const std::vector<gpcov> & covAllDimensionsInput,
                            const arma::vec & sigmaAllDimensionsInput,
                            const arma::vec & priorTemperatureInput,
                            const arma::mat & xInitInput,
                            const unsigned int numparamInput) {
    ThetaOptim objective(yobsInput, fOdeModelInput, covAllDimensionsInput, sigmaAllDimensionsInput, priorTemperatureInput, xInitInput, numparamInput);
    cppoptlib::LbfgsbSolver<ThetaOptim> solver;
    Eigen::VectorXd theta(numparamInput);
    theta.fill(1);
    solver.minimize(objective, theta);
    const arma::vec & thetaArgmin = arma::vec(theta.data(), numparamInput, false, false);
    return thetaArgmin;
}
