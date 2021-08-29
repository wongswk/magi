#define NDEBUG
#define BOOST_DISABLE_ASSERTS
#define EIGEN_NO_DEBUG

#include <armadillo>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgsbsolver.h>

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
    double freq = 1;
    const arma::vec & zmodEffectiveSq = arma::square(zmodEffective);
    freq = arma::sum(arma::linspace<arma::vec>(1, zmodEffective.size(), zmodEffective.size()) % zmodEffectiveSq) / arma::sum(zmodEffectiveSq);

    double meanFactor = 0.5 / freq;
    double sdFactor = (1 - meanFactor) / 3;

    return arma::vec({meanFactor, sdFactor});
}


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
arma::cube calcMeanCurve(const arma::vec & xInput,
                        const arma::vec & yInput,
                        const arma::vec & xOutput,
                        const arma::mat & phiCandidates,
                        const arma::vec & sigmaCandidates,
                        const std::string kerneltype = "generalMatern",
                        const bool useDeriv = false) {
    if(kerneltype != "generalMatern") std::cerr << "kerneltype other than generalMatern is not supported\n";

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
    const arma::vec & thetaArgmin = arma::vec(theta.data(), fOdeModelInput.thetaSize, true, false);
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
    const arma::mat & phiFull;
    const arma::uvec & missingComponentDim;

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
//        std::cout << "average priorFactor in PhiOptim =\n" << priorFactor << "\n";

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
    const arma::mat & phiArgmin = arma::mat(phi.data(), 2, missingComponentDim.size(), true, false);
    return phiArgmin;
}


class XmissingThetaPhiOptim : public cppoptlib::BoundedProblem<double> {
public:
    const arma::mat & yobs;
    const arma::vec & tvec;
    const OdeSystem & fOdeModel;
    const arma::vec & sigmaAllDimensions;
    const arma::vec & priorTemperature;
    arma::mat xInit;
    arma::vec thetaInit;
    arma::mat phiAllDimensions;
    const arma::uvec & missingComponentDim;
    const double SCALE = 1;

    double value(const Eigen::VectorXd & xthetaphiInput) override {
        if ((xthetaphiInput.array() < this->lowerBound().array()).any()){
            return INFINITY;
        }
        if ((xthetaphiInput.array() > this->upperBound().array()).any()){
            return INFINITY;
        }
        if (xthetaphiInput.array().isNaN().any()){
            return INFINITY;
        }

        const arma::vec xthetaphi(
                const_cast<double*>(xthetaphiInput.data()),
                xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * missingComponentDim.size(),
                false,
                false);

        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            xInit.col(missingComponentDim(id)) = xthetaphi.subvec(
                    xInit.n_rows * (id), xInit.n_rows * (id + 1) - 1);
        }

        thetaInit = xthetaphi.subvec(
                xInit.n_rows * missingComponentDim.size(), xInit.n_rows * missingComponentDim.size() + thetaInit.size() - 1);

        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            phiAllDimensions.col(missingComponentDim(id)) = xthetaphi.subvec(
                    xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * id,
                    xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * (id + 1) - 1);
//            std::cout << "value: phiAllDimensions.col(missingComponentDim(id)) = " << phiAllDimensions.col(missingComponentDim(id)) << "\n";
        }
//        std::cout << "inside evaluation of value\n";

        const lp & out = xthetaphisigmallik( xInit,
                                             thetaInit,
                                             phiAllDimensions,
                                             sigmaAllDimensions,
                                             yobs,
                                             tvec,
                                             fOdeModel);
        if (out.gradient.has_nan() || isnan(out.value)){
            return INFINITY;
        }
        return -out.value*SCALE;
    }

    void gradient(const Eigen::VectorXd & xthetaphiInput, Eigen::VectorXd & grad) override {
        if ((xthetaphiInput.array() < this->lowerBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < xthetaphiInput.size(); i++){
                if(xthetaphiInput[i] < this->lowerBound()[i]){
                    grad[i] = -1;
                }
            }
            return;
        }
        if ((xthetaphiInput.array() > this->upperBound().array()).any()){
            grad.fill(0);
            for(unsigned i = 0; i < xthetaphiInput.size(); i++){
                if(xthetaphiInput[i] < this->upperBound()[i]){
                    grad[i] = 1;
                }
            }
            return;
        }

        const arma::vec xthetaphi(
                const_cast<double*>(xthetaphiInput.data()),
                xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * missingComponentDim.size(),
                false,
                false);

        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            xInit.col(missingComponentDim(id)) = xthetaphi.subvec(
                    xInit.n_rows * (id), xInit.n_rows * (id + 1) - 1);
        }

        thetaInit = xthetaphi.subvec(
                xInit.n_rows * missingComponentDim.size(), xInit.n_rows * missingComponentDim.size() + thetaInit.size() - 1);

        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            phiAllDimensions.col(missingComponentDim(id)) = xthetaphi.subvec(
                    xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * id,
                    xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * (id + 1) - 1);
//            std::cout << "gradient: phiAllDimensions.col(missingComponentDim(id)) = " << phiAllDimensions.col(missingComponentDim(id)) << "\n";
        }

        lp out = xthetaphisigmallik( xInit,
                                             thetaInit,
                                             phiAllDimensions,
                                             sigmaAllDimensions,
                                             yobs,
                                             tvec,
                                             fOdeModel);
        out.gradient *= SCALE;
        out.value *= SCALE;

        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            for (unsigned j = 0; j < xInit.n_rows; j++){
                grad[xInit.n_rows * (id) + j] = -out.gradient(xInit.n_rows * missingComponentDim(id) + j);
            }
        }
        for (unsigned j = 0; j < thetaInit.size(); j++){
            grad[xInit.n_rows * missingComponentDim.size() + j] = -out.gradient(xInit.size() + j);
        }
        for (unsigned id = 0; id < missingComponentDim.size(); id++){
            for (unsigned j = 0; j < phiAllDimensions.n_rows; j++){
                grad[xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * id + j] =
                        -out.gradient(xInit.size() + thetaInit.size() + phiAllDimensions.n_rows * missingComponentDim(id) + j);
            }
        }
//        std::cout << "after gradient assignment =\n" << grad.transpose();
    }

    XmissingThetaPhiOptim(const arma::mat & yobsInput,
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
            phiAllDimensions(phiFullInput),
            missingComponentDim(missingComponentDimInput) {
        Eigen::VectorXd lb(xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * missingComponentDim.size());
        lb.fill(-INFINITY);
        Eigen::VectorXd ub(xInit.n_rows * missingComponentDim.size() + thetaInit.size() + phiAllDimensions.n_rows * missingComponentDim.size());
        ub.fill(INFINITY);

        for (unsigned j = 0; j < thetaInit.size(); j++){
            lb[xInit.n_rows * missingComponentDim.size() + j] = fOdeModel.thetaLowerBound(j) + 1e-6;
            ub[xInit.n_rows * missingComponentDim.size() + j] = fOdeModel.thetaUpperBound(j) - 1e-6;
        }

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
//        std::cout << "average priorFactor in PhiOptim =\n" << priorFactor << "\n";

        for(unsigned i = 0; i < missingComponentDim.size(); i++){
            ub[xInit.n_rows * missingComponentDim.size() + thetaInit.size() + 2*i] = maxScale * 5;
            lb[xInit.n_rows * missingComponentDim.size() + thetaInit.size() + 2*i] = maxScale * 1e-3;
            ub[xInit.n_rows * missingComponentDim.size() + thetaInit.size() + 2*i+1] = maxDist * 5;
            lb[xInit.n_rows * missingComponentDim.size() + thetaInit.size() + 2*i+1] = std::min(maxDist * priorFactor(0) * 0.5, minDist);
        }
        this->setLowerBound(lb);
        this->setUpperBound(ub);
//        std::cout << "finish set up of the problem\n";
//        std::cout << "ub = \n" << ub << "\nlb = \n" << lb << "\n";
    }
};

arma::mat optimizeXmissingThetaPhi(const arma::mat & yobsInput,
                                   const arma::vec & tvecInput,
                                   const OdeSystem & fOdeModelInput,
                                   const arma::vec & sigmaAllDimensionsInput,
                                   const arma::vec & priorTemperatureInput,
                                   const arma::mat & xInitInput,
                                   const arma::vec & thetaInitInput,
                                   const arma::mat & phiInitInput,
                                   const arma::uvec & missingComponentDim) {
    XmissingThetaPhiOptim objective(yobsInput, tvecInput, fOdeModelInput, sigmaAllDimensionsInput, priorTemperatureInput, xInitInput, thetaInitInput, phiInitInput, missingComponentDim);
    cppoptlib::LbfgsbSolver<XmissingThetaPhiOptim> solver;

    Eigen::VectorXd xThetaPhi(xInitInput.n_rows * missingComponentDim.size() + thetaInitInput.size() + phiInitInput.n_rows * missingComponentDim.size());

    for (unsigned id = 0; id < missingComponentDim.size(); id++){
        for (unsigned j = 0; j < xInitInput.n_rows; j++){
            xThetaPhi[xInitInput.n_rows * (id) + j] = xInitInput(j, missingComponentDim(id));
        }
    }
    for (unsigned j = 0; j < thetaInitInput.size(); j++){
        xThetaPhi[xInitInput.n_rows * missingComponentDim.size() + j] = thetaInitInput(j);
    }
    for (unsigned id = 0; id < missingComponentDim.size(); id++){
        for (unsigned j = 0; j < phiInitInput.n_rows; j++){
            xThetaPhi[xInitInput.n_rows * missingComponentDim.size() + thetaInitInput.size() + phiInitInput.n_rows * id + j] =
                    phiInitInput(j, missingComponentDim(id));
        }
    }

//    std::cout << "inside optimizeXmissingThetaPhi\n"
//    << "init xThetaPhi = " << xThetaPhi;
    Eigen::VectorXd xThetaPhiInit = xThetaPhi;

    solver.minimize(objective, xThetaPhi);
    if (! isfinite(objective.value(xThetaPhi))){
        const arma::vec & xThetaPhiArgmin = arma::vec(xThetaPhiInit.data(), xInitInput.n_rows * missingComponentDim.size() + thetaInitInput.size() + phiInitInput.n_rows * missingComponentDim.size(), true, false);
        return xThetaPhiArgmin;
    }

    if (objective.value(xThetaPhi) < objective.value(xThetaPhiInit)){
        const arma::vec & xThetaPhiArgmin = arma::vec(xThetaPhi.data(), xInitInput.n_rows * missingComponentDim.size() + thetaInitInput.size() + phiInitInput.n_rows * missingComponentDim.size(), true, false);
        return xThetaPhiArgmin;
    }else{
        const arma::vec & xThetaPhiArgmin = arma::vec(xThetaPhiInit.data(), xInitInput.n_rows * missingComponentDim.size() + thetaInitInput.size() + phiInitInput.n_rows * missingComponentDim.size(), true, false);
        return xThetaPhiArgmin;
    }
}
