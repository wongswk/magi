#include <armadillo>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgsbsolver.h>
#include <cppoptlib/solver/bfgssolver.h>

#include "tgtdistr.h"
#include "fullloglikelihood.h"


// [[Rcpp::export]]
arma::vec calcFrequencyBasedPrior(const arma::vec & x){
    arma::cx_mat z = arma::fft(x);
    arma::vec zmod(z.n_rows);
    for(unsigned int i = 0; i < z.n_rows; i++){
        zmod(i) = std::sqrt(std::norm(z(i, 0)));
    }
    const arma::vec & zmodEffective = zmod(arma::span(1, (zmod.size() - 1) / 2));
    const arma::vec zmodEffectiveSorted = arma::sort(zmodEffective);
    double upperQuarter = zmodEffectiveSorted(static_cast<unsigned int>(std::ceil(zmodEffectiveSorted.size() * 0.75)) - 1);
    double lowerQuarter = zmodEffectiveSorted(std::max(static_cast<int>(std::floor(zmodEffectiveSorted.size() * 0.25)) - 1, 0));
    double iqr = upperQuarter - lowerQuarter;
    arma::vec outliers = zmodEffective(arma::find(zmodEffective > upperQuarter + 1.5 * iqr));
    long long int freq = 1;
    if(!outliers.empty()){
        for(unsigned long long int i = zmodEffective.size(); i > 0; i--){
            if(arma::min(arma::abs(zmodEffective(i - 1) - outliers)) < 1e-6){
                freq = i;
                break;
            }
        }

    }else{
        freq = zmodEffective.index_max() + 1;
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
    const double sigmaExogenScalar;
    const bool useFrequencyBasedPrior;
    arma::vec priorFactor;
    double maxDist;

    double value(const Eigen::VectorXd & phisigInput) override {
        if ((phisigInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        if ((phisigInput.array() > this->upperBound().array()).any()){
            return INFINITY;
        }
        arma::vec phisig = arma::vec(const_cast<double*>(phisigInput.data()), numparam, false, false);
        if(sigmaExogenScalar > 0){
            phisig = arma::join_vert(phisig, arma::vec({sigmaExogenScalar}));
        }
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
        if ((phisigInput.array() > this->upperBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < numparam; i++){
                if(phisigInput[i] > this->upperBound()[i]){
                    grad[i] = 1;
                }
            }
            return;
        }
        arma::vec phisig = arma::vec(const_cast<double*>(phisigInput.data()), numparam, false, false);
        if(sigmaExogenScalar > 0){
            phisig = arma::join_vert(phisig, arma::vec({sigmaExogenScalar}));
        }
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
                                const double sigmaExogenScalarInput,
                                const bool useFrequencyBasedPriorInput) :
            BoundedProblem(numparamInput),
            kernel(std::move(kernelInput)),
            yobs(yobsInput),
            dist(distInput),
            numparam(numparamInput),
            sigmaExogenScalar(sigmaExogenScalarInput),
            useFrequencyBasedPrior(useFrequencyBasedPriorInput) {
        std::cout << "kernelInput in PhiGaussianProcessSmoothing is " << kernelInput << "\n";
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

        Eigen::VectorXd lb(numparam);
        lb.fill(1e-4);
        this->setLowerBound(lb);

        maxDist = dist.max();
        double maxScale = arma::max(arma::abs(yobs(arma::find_finite(yobs))));
        maxScale = std::max(maxScale, maxDist);

        Eigen::VectorXd ub(numparam);
        ub.fill(10 * maxScale);
        for(unsigned i = 0; i < yobsInput.n_cols; i++) {
            const arma::uvec finite_elem = arma::find_finite(yobs.col(i));
            if (finite_elem.size() > 0){
                ub[phiDim * i] = arma::max(arma::abs((yobs.col(i).eval().elem(finite_elem))));
            }
            ub[phiDim * i + 1] = maxDist;
        }
        this->setUpperBound(ub);

        priorFactor = arma::zeros(2);
        if(useFrequencyBasedPrior){
            for (unsigned j = 0; j < yobs.n_cols; j++){
                priorFactor += calcFrequencyBasedPrior(yobs.col(j));
            }
            priorFactor /= yobs.n_cols;
            std::cout << "priorFactor =\n" << priorFactor << "\n";
        }
    }
};


// [[Rcpp::export]]
arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   const double sigmaExogenScalar = -1,
                   bool useFrequencyBasedPrior = false) {
    std::cout << "kernelInput in gpsmooth is " << kernelInput << "\n";
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

    PhiGaussianProcessSmoothing objective(yobsInput, distInput, std::move(kernelInput), numparam, sigmaExogenScalar, useFrequencyBasedPrior);
    cppoptlib::LbfgsbSolver<PhiGaussianProcessSmoothing> solver;
    Eigen::VectorXd phisig(numparam);
    phisig.fill(1);
    double maxDist = distInput.max();
    double sdOverall = 0;
    for(unsigned i = 0; i < yobsInput.n_cols; i++) {
        phisig[phiDim * i] = 0.5 * arma::stddev(yobsInput.col(i));
        phisig[phiDim * i + 1] = 0.5 * maxDist;
        sdOverall += phisig[phiDim * i];
    }
    if(sigmaExogenScalar <= 0){
        phisig[phiDim * yobsInput.n_cols] = sdOverall / yobsInput.n_cols;
    }
    solver.minimize(objective, phisig);
    const arma::vec & phisigArgmin = arma::vec(phisig.data(), numparam, false, false);
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
    const bool useBand;

    double value(const Eigen::VectorXd & thetaInput) override {
        if ((thetaInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        if ((thetaInput.array() > this->upperBound().array()).any()){
            return INFINITY;
        }
        const arma::vec & xtheta = arma::join_vert(
            arma::vectorise(xInit),
            arma::vec(const_cast<double*>(thetaInput.data()), fOdeModel.thetaSize, false, false)
        );
        const lp & out = xthetallik(
                xtheta,
                covAllDimensions,
                sigmaAllDimensions,
                yobs,
                fOdeModel,
                useBand,
                priorTemperature);
        return -out.value;
    }

    void gradient(const Eigen::VectorXd & thetaInput, Eigen::VectorXd & grad) override {
        if ((thetaInput.array() < this->lowerBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < fOdeModel.thetaSize; i++){
                if(thetaInput[i] < this->lowerBound()[i]){
                    grad[i] = -1;
                }
            }
            return;
        }
        if ((thetaInput.array() > this->upperBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < fOdeModel.thetaSize; i++){
                if(thetaInput[i] > this->upperBound()[i]){
                    grad[i] = 1;
                }
            }
            return;
        }
        const arma::vec & xtheta = arma::join_vert(
                arma::vectorise(xInit),
                arma::vec(const_cast<double*>(thetaInput.data()), fOdeModel.thetaSize, false, false)
        );
        const lp & out = xthetallik(
                xtheta,
                covAllDimensions,
                sigmaAllDimensions,
                yobs,
                fOdeModel,
                useBand,
                priorTemperature);
        for(unsigned i = 0; i < fOdeModel.thetaSize; i++){
            grad[i] = -out.gradient(xInit.size() + i);
        }
    }

    ThetaOptim(const arma::mat & yobsInput,
               const OdeSystem & fOdeModelInput,
               const std::vector<gpcov> & covAllDimensionsInput,
               const arma::vec & sigmaAllDimensionsInput,
               const arma::vec & priorTemperatureInput,
               const arma::mat & xInitInput,
               const bool useBandInput) :
            BoundedProblem(fOdeModelInput.thetaSize),
            yobs(yobsInput),
            fOdeModel(fOdeModelInput),
            covAllDimensions(covAllDimensionsInput),
            sigmaAllDimensions(sigmaAllDimensionsInput),
            priorTemperature(priorTemperatureInput),
            xInit(xInitInput),
            useBand(useBandInput) {
        const Eigen::Map<Eigen::VectorXd> lb (const_cast<double*>(fOdeModel.thetaLowerBound.memptr()), fOdeModel.thetaSize);
        this->setLowerBound(lb.array() + 1e-6);
        const Eigen::Map<Eigen::VectorXd> ub (const_cast<double*>(fOdeModel.thetaUpperBound.memptr()), fOdeModel.thetaSize);
        this->setUpperBound(ub.array() - 1e-6);
    }
};

// [[Rcpp::export]]
arma::vec optimizeThetaInit(const arma::mat & yobsInput,
                            const OdeSystem & fOdeModelInput,
                            const std::vector<gpcov> & covAllDimensionsInput,
                            const arma::vec & sigmaAllDimensionsInput,
                            const arma::vec & priorTemperatureInput,
                            const arma::mat & xInitInput,
                            const bool useBandInput) {
    ThetaOptim objective(yobsInput, fOdeModelInput, covAllDimensionsInput, sigmaAllDimensionsInput, priorTemperatureInput, xInitInput, useBandInput);
    cppoptlib::LbfgsbSolver<ThetaOptim> solver;
    Eigen::VectorXd theta(fOdeModelInput.thetaSize);
    theta.fill(1);
    solver.minimize(objective, theta);
    const arma::vec & thetaArgmin = arma::vec(theta.data(), fOdeModelInput.thetaSize, false, false);
    return thetaArgmin;
}


class PhiOptim : public cppoptlib::BoundedProblem<double> {
public:
    const arma::mat & yobs;
    const arma::vec & tvec;
    const OdeSystem & fOdeModel;
    const arma::vec & sigmaAllDimensions;
    const arma::vec & priorTemperature;
    const arma::mat & xInit;
    const arma::vec & thetaInit;
    const arma::uvec & missingComponentDim;
    const arma::mat & phiFull;

    double value(const Eigen::VectorXd & phiInput) override {
        if ((phiInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        const arma::mat phiMissingDimensions(
                const_cast<double*>(phiInput.data()),
                2,
                missingComponentDim.size(),
                false,
                false);
        arma::mat phiAllDimensions = phiFull;
        phiAllDimensions.cols(missingComponentDim) = phiMissingDimensions;
        const lp & out = xthetaphisigmallik( xInit,
                                             thetaInit,
                                             phiAllDimensions,
                                             sigmaAllDimensions,
                                             yobs,
                                             tvec,
                                             fOdeModel);
        return -out.value;
    }

    void gradient(const Eigen::VectorXd & phiInput, Eigen::VectorXd & grad) override {
        if ((phiInput.array() < this->lowerBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < phiInput.size(); i++){
                if(phiInput[i] < this->lowerBound()[i]){
                    grad[i] = -1;
                }
            }
            return;
        }
        const arma::mat phiMissingDimensions(
                const_cast<double*>(phiInput.data()),
                2,
                missingComponentDim.size(),
                false,
                false);
        arma::mat phiAllDimensions = phiFull;
        phiAllDimensions.cols(missingComponentDim) = phiMissingDimensions;

        const lp & out = xthetaphisigmallik( xInit,
                                             thetaInit,
                                             phiAllDimensions,
                                             sigmaAllDimensions,
                                             yobs,
                                             tvec,
                                             fOdeModel);

        for(unsigned i = 0; i < missingComponentDim.size(); i++){
            unsigned currentDim = missingComponentDim[i];
            grad[2*i] = -out.gradient(xInit.size() + thetaInit.size() + 2*currentDim);
            grad[2*i+1] = -out.gradient(xInit.size() + thetaInit.size() + 2*currentDim + 1);
        }
    }

    PhiOptim(const arma::mat & yobsInput,
             const arma::vec & tvecInput,
             const OdeSystem & fOdeModelInput,
             const arma::vec & sigmaAllDimensionsInput,
             const arma::vec & priorTemperatureInput,
             const arma::mat & xInitInput,
             const arma::vec & thetaInitInput,
             const arma::mat & phiFullInput,
             const arma::uvec & missingComponentDimInput) :
            BoundedProblem(missingComponentDimInput.size() * 2),
            yobs(yobsInput),
            tvec(tvecInput),
            fOdeModel(fOdeModelInput),
            sigmaAllDimensions(sigmaAllDimensionsInput),
            priorTemperature(priorTemperatureInput),
            xInit(xInitInput),
            thetaInit(thetaInitInput),
            phiFull(phiFullInput),
            missingComponentDim(missingComponentDimInput) {
        Eigen::VectorXd lb(missingComponentDim.size() * 2);
        Eigen::VectorXd ub(missingComponentDim.size() * 2);

        const double maxDist = (arma::max(tvecInput) - arma::min(tvecInput));
        const double minDist = arma::min(arma::abs(arma::diff(tvecInput)));
        const double maxScale = arma::max(arma::abs(yobs(arma::find_finite(yobs))));

        arma::vec priorFactor = arma::zeros(2);
        for (unsigned j = 0; j < yobs.n_cols; j++){
            if (arma::any(missingComponentDim == j)){
                continue;
            }
            const arma::vec & yobsThisDim = yobs.col(j);
            priorFactor += calcFrequencyBasedPrior(yobsThisDim(arma::find_finite(yobsThisDim)));
        }
        priorFactor /= (yobs.n_cols - missingComponentDim.size());
        std::cout << "average priorFactor in PhiOptim =\n" << priorFactor << "\n";

        for(unsigned i = 0; i < missingComponentDim.size(); i++){
            ub[2*i] = maxScale * 5;
            lb[2*i] = maxScale * 1e-3;
            ub[2*i+1] = maxDist * 5;
            lb[2*i+1] = std::min(maxDist * priorFactor(0) * 0.5, minDist);
        }
        this->setLowerBound(lb);
        this->setUpperBound(ub);
    }
};


// [[Rcpp::export]]
arma::mat optimizePhi(const arma::mat & yobsInput,
                      const arma::vec & tvecInput,
                      const OdeSystem & fOdeModelInput,
                      const arma::vec & sigmaAllDimensionsInput,
                      const arma::vec & priorTemperatureInput,
                      const arma::mat & xInitInput,
                      const arma::vec & thetaInitInput,
                      const arma::mat & phiInitInput,
                      const arma::uvec & missingComponentDim) {
    PhiOptim objective(yobsInput, tvecInput, fOdeModelInput, sigmaAllDimensionsInput, priorTemperatureInput, xInitInput, thetaInitInput, phiInitInput, missingComponentDim);
    cppoptlib::LbfgsbSolver<PhiOptim> solver;
    Eigen::VectorXd phi(2 * missingComponentDim.size());
    for(unsigned i = 0; i < missingComponentDim.size(); i++){
        unsigned currentDim = missingComponentDim[i];
        phi[2*i] = phiInitInput(0, currentDim);
        phi[2*i+1] = phiInitInput(1, currentDim);
    }
    solver.minimize(objective, phi);
    const arma::mat & phiArgmin = arma::mat(phi.data(), 2, missingComponentDim.size(), false, false);
    return phiArgmin;
}
