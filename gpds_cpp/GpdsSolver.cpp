#include "GpdsSolver.h"
#include "gpsmoothing.h"
#include "tgtdistr.h"
#include "fullloglikelihood.h"
#include "Sampler.h"


GpdsSolver::GpdsSolver(const arma::mat & yFull,
                       const OdeSystem & odeModel,
                       const arma::vec & tvecFull,
                       const arma::vec & sigmaExogenous,
                       const double priorTemperatureLevel,
                       const double priorTemperatureDeriv,
                       std::string kernel,
                       const int nstepsHmc,
                       const double burninRatioHmc,
                       const unsigned int niterHmc,
                       const double stepSizeFactorHmc,
                       const int nEpoch,
                       const int bandSize,
                       bool useFrequencyBasedPrior,
                       bool useBand,
                       bool useMean,
                       bool useScalerSigma,
                       bool useFixedSigma,
                       bool verbose) :
        yFull(yFull),
        odeModel(odeModel),
        tvecFull(tvecFull),
        sigmaExogenous(sigmaExogenous),
        priorTemperature({priorTemperatureLevel, priorTemperatureDeriv}),
        kernel(kernel),
        nstepsHmc(nstepsHmc),
        burninRatioHmc(burninRatioHmc),
        niterHmc(niterHmc),
        stepSizeFactorHmc(stepSizeFactorHmc),
        nEpoch(nEpoch),
        bandSize(bandSize),
        useFrequencyBasedPrior(useFrequencyBasedPrior),
        useBand(useBand),
        useMean(useMean),
        useScalerSigma(useScalerSigma),
        useFixedSigma(useFixedSigma),
        verbose(verbose),
        ydim(yFull.n_cols),
        sigmaSize(useScalerSigma ? 1 : yFull.n_cols),
        distSignedFull(tvecFull.size(), tvecFull.size()),
        indicatorRowWithObs(yFull.n_rows),
        indicatorMatWithObs(yFull.n_rows, yFull.n_cols, arma::fill::zeros),
        phiAllDimensions(2, yFull.n_cols)
{
    if(kernel != "generalMatern"){
        throw std::runtime_error("only generalMatern kernel has full support");
    }
    // generate intermediate data
    if(useBand && useMean){
        loglikflag = "withmeanBand";
    }else if(useBand && !useMean){
        loglikflag = "band";
    }else if(!useBand && useMean){
        loglikflag = "withmean";
    }else if(!useBand && !useMean){
        loglikflag = "usual";
    }

    for(unsigned int i = 0; i < distSignedFull.n_cols; i++){
        distSignedFull.col(i) = tvecFull - tvecFull(i);
    }


    for(unsigned int i = 0; i < yFull.n_rows; i++) {
        for(unsigned int j = 0; j < ydim; j++) {
            indicatorMatWithObs(i, j) = std::isfinite(yFull(i, j));
        }
    }
    indicatorRowWithObs = arma::sum(indicatorMatWithObs, 1);
    idxRowWithObs = arma::find(indicatorRowWithObs > 0);
    yObs = yFull.rows(idxRowWithObs);
    distSignedObs = distSignedFull.submat(idxRowWithObs, idxRowWithObs);
    idxColElemWithObs.resize(ydim);
    for(unsigned int j = 0; j < ydim; j++) {
        idxColElemWithObs[j] = arma::find(indicatorMatWithObs.col(j) > 0);
    }

    covAllDimensions.resize(ydim);

    if(kernel == "matern"){
        kernelCov = maternCov;
    }else if(kernel == "rbf"){
        kernelCov = rbfCov;
    }else if(kernel == "compact1"){
        kernelCov = compact1Cov;
    }else if(kernel == "periodicMatern"){
        kernelCov = periodicMaternCov;
    }else if(kernel == "generalMatern"){
        kernelCov = generalMaternCov;
    }else{
        throw std::runtime_error("kernel is not specified correctly");
    }

}

void GpdsSolver::setupPhiSigma() {
    if(sigmaExogenous.empty()){
        if(useScalerSigma){
            const arma::vec & phisig = gpsmooth(yObs,
                                                arma::abs(distSignedObs),
                                                kernel,
                                                -1,
                                                useFrequencyBasedPrior);
            phiAllDimensions = arma::reshape(phisig.subvec(0, 2*ydim - 1).eval(), 2, ydim);
            sigmaInit = phisig.subvec(2*ydim, 2*ydim);
        }else{
            sigmaInit.resize(ydim);
            arma::uvec sucess(ydim);
            for(unsigned j = 0; j < ydim; j++){
                if(idxColElemWithObs[j].size() >= 3){
                    const arma::vec & yObsCol = yFull.col(j).eval().elem(idxColElemWithObs[j]);
                    const arma::mat & distSignedObsCol = distSignedFull.submat(idxColElemWithObs[j], idxColElemWithObs[j]);
                    const arma::vec & phisig = gpsmooth(yObsCol,
                                                        arma::abs(distSignedObsCol),
                                                        kernel,
                                                        -1,
                                                        useFrequencyBasedPrior);
                    phiAllDimensions.col(j) = phisig.subvec(0, 1);
                    sigmaInit(j) = phisig(2);
                    sucess(j) = 1;
                }else{
                    sucess(j) = 0;
                }
            }
            const arma::uvec & sucessDim = arma::find(sucess > 0);
            const arma::vec & phiMean = arma::mean(phiAllDimensions.cols(sucessDim), 1);
            const double sigmaMean = arma::mean(sigmaInit(sucessDim));

            const arma::uvec & failedDim = arma::find(sucess == 0);
            sigmaInit(failedDim).fill(sigmaMean);
            phiAllDimensions.submat(arma::uvec({0}), failedDim).fill(phiMean(0));
            phiAllDimensions.submat(arma::uvec({1}), failedDim).fill(phiMean(1));
        }
    }else{
        if(useScalerSigma){
            const arma::vec & phisig = gpsmooth(yObs,
                                                arma::abs(distSignedObs),
                                                kernel,
                                                sigmaExogenous(0),
                                                useFrequencyBasedPrior);
            phiAllDimensions = arma::reshape(phisig.subvec(0, 2*ydim - 1).eval(), 2, ydim);
            sigmaInit = sigmaExogenous.subvec(0, 0);
        }else{
            sigmaInit = sigmaExogenous;
            arma::uvec sucess(ydim);
            for(unsigned j = 0; j < ydim; j++){
                if(idxColElemWithObs[j].size() >= 3){
                    const arma::vec & yObsCol = yFull.col(j).eval().elem(idxColElemWithObs[j]);
                    const arma::mat & distSignedObsCol = distSignedFull.submat(idxColElemWithObs[j], idxColElemWithObs[j]);
                    const arma::vec & phisig = gpsmooth(yObsCol,
                                                        arma::abs(distSignedObsCol),
                                                        kernel,
                                                        sigmaExogenous(j),
                                                        useFrequencyBasedPrior);
                    phiAllDimensions.col(j) = phisig.subvec(0, 1);
                    sucess(j) = 1;
                }else{
                    sucess(j) = 0;
                }
            }
            const arma::uvec & sucessDim = arma::find(sucess > 0);
            const arma::vec & phiMean = arma::mean(phiAllDimensions.cols(sucessDim), 1);

            const arma::uvec & failedDim = arma::find(sucess == 0);
            phiAllDimensions.submat(arma::uvec({0}), failedDim).fill(phiMean(0));
            phiAllDimensions.submat(arma::uvec({1}), failedDim).fill(phiMean(1));
        }
    }

    for(unsigned j = 0; j < ydim; j++){
        covAllDimensions[j] = kernelCov(phiAllDimensions.col(j), distSignedFull, 3);
        covAllDimensions[j].addBandCov(bandSize);
    }
}

void GpdsSolver::initXmudotmu() {
    arma::vec sigmaUsed(ydim);
    if(useScalerSigma){
        sigmaUsed.fill(sigmaInit(0));
    }else{
        sigmaUsed = sigmaInit;
    }
    arma::uvec sucess(ydim);
    arma::mat xInitAllDim = arma::ones(yFull.n_rows, ydim);
    for(unsigned j = 0; j < ydim; j++){
        if(idxColElemWithObs[j].size() >= 3){
            const arma::vec & yObsCol = yFull.col(j).eval().elem(idxColElemWithObs[j]);
            const arma::vec & tvecObsCol = tvecFull.elem(idxColElemWithObs[j]);

            const arma::cube & xdx = calcMeanCurve(tvecObsCol,
                                                   yObsCol,
                                                   tvecFull,
                                                   phiAllDimensions.col(j),
                                                   sigmaUsed.subvec(j, j),
                                                   kernel,
                                                   true);
            xInitAllDim.col(j) = xdx.slice(0);

            covAllDimensions[j].mu = xdx.slice(0);
            covAllDimensions[j].dotmu = xdx.slice(1);
            sucess(j) = 2;
        }else if(!idxColElemWithObs[j].empty()){
            const arma::vec & yObsCol = yFull.col(j).eval().elem(idxColElemWithObs[j]);
            xInitAllDim.col(j).fill(arma::mean(yObsCol));
//            xInitAllDim.col(j) = arma::randn(tvecFull.size()) + arma::mean(yObsCol);

            covAllDimensions[j].mu = arma::ones(tvecFull.size()) * arma::mean(yObsCol);
            covAllDimensions[j].dotmu = arma::zeros(tvecFull.size());
            sucess(j) = 1;
        }else{
            xInitAllDim.col(j) = arma::ones(tvecFull.size());
//            xInitAllDim.col(j) = arma::randn(yFull.n_rows) + 1;

            covAllDimensions[j].mu = arma::zeros(tvecFull.size());
            covAllDimensions[j].dotmu = arma::zeros(tvecFull.size());
            sucess(j) = 0;
        }
    }

    const arma::uvec & sucessDim = arma::find(sucess == 2);
    const arma::vec & xMean = arma::mean(xInitAllDim.cols(sucessDim), 1);
    const arma::uvec & failedDim1 = arma::find(sucess == 1);

    xInit = xInitAllDim;
}

void GpdsSolver::initTheta() {
    arma::vec sigmaUsed(ydim);
    if(useScalerSigma){
        sigmaUsed.fill(sigmaInit(0));
    }else{
        sigmaUsed = sigmaInit;
    }

    thetaInit = optimizeThetaInit(yFull,
                                  odeModel,
                                  covAllDimensions,
                                  sigmaUsed,
                                  priorTemperature,
                                  xInit);
}

void GpdsSolver::initMissingComponent() {
    const unsigned int nSGD = 100;
    double learningRate = 1e-8;
    const arma::uvec & nobsEachDim = arma::sum(indicatorMatWithObs, 0).t();
    const arma::uvec & missingComponentDim = arma::find(nobsEachDim < 3);
    if(missingComponentDim.empty()){
        return;
    }

    const arma::uvec & observedComponentDim = arma::find(nobsEachDim >= 0);
    for (auto iPtr = missingComponentDim.begin(); iPtr < missingComponentDim.end(); iPtr++){
        xInit.col(*iPtr) = arma::mean(xInit.cols(observedComponentDim), 1);
        xInit.submat(idxColElemWithObs[*iPtr], arma::uvec({*iPtr})) =
                yFull.submat(idxColElemWithObs[*iPtr], arma::uvec({*iPtr}));
    }

    for(unsigned int iSGD = 0; iSGD < nSGD; iSGD++){
        const lp & llik = xthetaphisigmallik( xInit,
                                              thetaInit,
                                              phiAllDimensions,
                                              sigmaInit,
                                              yFull,
                                              tvecFull,
                                              odeModel);
        thetaInit += learningRate * llik.gradient.subvec(xInit.size(), xInit.size() + thetaInit.size() - 1);
        for (auto iPtr = missingComponentDim.begin(); iPtr < missingComponentDim.end(); iPtr++){
            xInit.col(*iPtr) += learningRate * llik.gradient.subvec(
                    xInit.n_rows * (*iPtr), xInit.n_rows * (*iPtr + 1) - 1);
            phiAllDimensions.col(*iPtr) += learningRate * llik.gradient.subvec(
                    xInit.size() + thetaInit.size() + phiAllDimensions.n_rows * (*iPtr),
                    xInit.size() + thetaInit.size() + phiAllDimensions.n_rows * (*iPtr + 1) - 1);
        }
        thetaInit = arma::max(thetaInit, odeModel.thetaLowerBound + 1e-6);
        thetaInit = arma::min(thetaInit, odeModel.thetaUpperBound + 1e-6);
        phiAllDimensions = arma::max(phiAllDimensions, arma::ones(arma::size(phiAllDimensions)) * 1e-6);
        if (verbose) {
            std::cout << "initMissingComponent iteration " << iSGD
                      << "; xthetaphisigmallik = " << llik.value
                      << "; phi missing dim = \n" << phiAllDimensions.cols(missingComponentDim)
                      << "\n";
        }
    }
}

void GpdsSolver::doHMC() {
    Sampler hmcSampler(yFull,
                       covAllDimensions,
                       nstepsHmc,
                       loglikflag,
                       priorTemperature,
                       sigmaSize,
                       odeModel,
                       niterHmc,
                       burninRatioHmc);
    arma::vec xthetasigmaInit = arma::join_vert(arma::join_vert(arma::vectorise(xInit), thetaInit), sigmaInit);
    arma::vec stepLowInit(xthetasigmaInit.size());
    stepLowInit.fill(1.0 / nstepsHmc * stepSizeFactorHmc);
    if(!sigmaExogenous.empty()){
        stepLowInit.subvec(xInit.size() + thetaInit.size(), stepLowInit.size() - 1).fill(0);
    }
    hmcSampler.sampleChian(xthetasigmaInit, stepLowInit, verbose);
    llikSamples = hmcSampler.lliklist;
    xthetasigmaSamples = hmcSampler.xth;
    stepLow = hmcSampler.stepLow;
}
