#include "MagiSolver.h"
#include "gpsmoothing.h"
#include "tgtdistr.h"
#include "fullloglikelihood.h"
#include "Sampler.h"


MagiSolver::MagiSolver(const arma::mat & yFull,
                       const OdeSystem & odeModel,
                       const arma::vec & tvecFull,
                       const arma::vec & sigmaExogenous,
                       const arma::mat & phiExogenous,
                       const arma::mat & xInitExogenous,
                       const arma::vec & thetaInitExogenous,
                       const arma::mat & muExogenous,
                       const arma::mat & dotmuExogenous,
                       const double priorTemperatureLevel,
                       const double priorTemperatureDeriv,
                       const double priorTemperatureObs,
                       std::string kernel,
                       const int nstepsHmc,
                       const double burninRatioHmc,
                       const unsigned int niterHmc,
                       const arma::vec stepSizeFactorHmcInput,
                       const arma::cube distSignedFullInput,
                       const int nEpoch,
                       const int bandSize,
                       bool useFrequencyBasedPrior,
                       bool useBand,
                       bool useMean,
                       bool useScalerSigma,
                       bool useFixedSigma,
                       bool skipMissingComponentOptimization,
                       bool verbose) :
        yFull(yFull),
        odeModel(odeModel),
        tvecFull(tvecFull),
        sigmaExogenous(sigmaExogenous),
        phiExogenous(phiExogenous),
        xInitExogenous(xInitExogenous),
        thetaInitExogenous(thetaInitExogenous),
        muExogenous(muExogenous),
        dotmuExogenous(dotmuExogenous),
        priorTemperature({priorTemperatureLevel, priorTemperatureDeriv, priorTemperatureObs}),
        kernel(kernel),
        nstepsHmc(nstepsHmc),
        burninRatioHmc(burninRatioHmc),
        niterHmc(niterHmc),
        nEpoch(nEpoch),
        bandSize(bandSize),
        useFrequencyBasedPrior(useFrequencyBasedPrior),
        useBand(useBand),
        useMean(useMean),
        useScalerSigma(useScalerSigma),
        useFixedSigma(useFixedSigma),
        skipMissingComponentOptimization(skipMissingComponentOptimization),
        verbose(verbose),
        ydim(yFull.n_cols),
        sigmaSize(useScalerSigma ? 1 : yFull.n_cols),
        distSignedFull(tvecFull.size(), tvecFull.size(), yFull.n_cols),
        indicatorRowWithObs(yFull.n_rows),
        indicatorMatWithObs(yFull.n_rows, yFull.n_cols, arma::fill::zeros),
        phiAllDimensions(2, yFull.n_cols),
        llikxthetasigmaSamples(1 + yFull.size() + odeModel.thetaSize + sigmaSize, niterHmc, nEpoch)
{
    if (stepSizeFactorHmcInput.n_elem == 0){
        stepSizeFactorHmc = arma::ones(1);
    }else{
        stepSizeFactorHmc = stepSizeFactorHmcInput;
    }

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

    if(distSignedFullInput.empty()){
        for(unsigned int j = 0; j < ydim; j++){
            for(unsigned int i = 0; i < distSignedFull.n_cols; i++){
                distSignedFull.slice(j).col(i) = tvecFull - tvecFull(i);
            }
        }
    }else{
        distSignedFull = distSignedFullInput;
    }


    for(unsigned int i = 0; i < yFull.n_rows; i++) {
        for(unsigned int j = 0; j < ydim; j++) {
            indicatorMatWithObs(i, j) = std::isfinite(yFull(i, j));
        }
    }
    indicatorRowWithObs = arma::sum(indicatorMatWithObs, 1);
    idxRowWithObs = arma::find(indicatorRowWithObs > 0);
    yObs = yFull.rows(idxRowWithObs);
    distSignedObs = distSignedFull.slice(0).submat(idxRowWithObs, idxRowWithObs);
    for(unsigned int j = 1; j < ydim; j++) {
        distSignedObs += distSignedFull.slice(j).submat(idxRowWithObs, idxRowWithObs);
    }
    distSignedObs /= ydim;

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

void MagiSolver::setupPhiSigma() {
    if(sigmaExogenous.empty() && phiExogenous.empty()){
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
                    const arma::mat & distSignedObsCol = distSignedFull.slice(j).submat(idxColElemWithObs[j], idxColElemWithObs[j]);
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
    }else if(!sigmaExogenous.empty() && phiExogenous.empty()){
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
                    const arma::mat & distSignedObsCol = distSignedFull.slice(j).submat(idxColElemWithObs[j], idxColElemWithObs[j]);
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
    }else if(!sigmaExogenous.empty() && !phiExogenous.empty()) {
        if(useScalerSigma){
            sigmaInit = sigmaExogenous.subvec(0, 0);
        }else{
            sigmaInit = sigmaExogenous;
        }
        phiAllDimensions = phiExogenous;
    }else{
        throw std::runtime_error("when supplying phiExogenous, sigmaExogenous must be supplied");
    }

    for(unsigned j = 0; j < ydim; j++){
        covAllDimensions[j] = kernelCov(phiAllDimensions.col(j), distSignedFull.slice(j), 3);
        covAllDimensions[j].tvecCovInput = tvecFull;

        // Workaround for phi1 getting too large and matrix inverse failing
        while (covAllDimensions[j].Cinv.n_rows == 0 || covAllDimensions[j].Kinv.n_rows == 0) {
          std::cout << "Cinv or Kinv failed for component " << j << " with phi1 = " << phiAllDimensions(0,j) << " and phi2 = " << phiAllDimensions(1,j) << endl;
          phiAllDimensions(0,j) *= 0.8;
          covAllDimensions[j] = kernelCov(phiAllDimensions.col(j), distSignedFull.slice(j), 3);
        }
        
        covAllDimensions[j].addBandCov(bandSize);
        
        // Diagnostic information
//        std::cout << "Component " << j << " Cinv max element: " << arma::max(arma::max(arma::abs(covAllDimensions[j].Cinv))) << ", Kinv max element: " << arma::max(arma::max(arma::abs(covAllDimensions[j].Kinv))) << endl;
//        std::cout << "Component " << j << " Cinv min element: " << arma::min(arma::min(arma::abs(covAllDimensions[j].Cinv))) << ", Kinv min element: " << arma::min(arma::min(arma::abs(covAllDimensions[j].Kinv))) << endl;

        
    }
}

void MagiSolver::initXmudotmu() {
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

    xInit = xInitAllDim;
    if(!xInitExogenous.empty()){
        xInit = xInitExogenous;
    }
    if(!muExogenous.empty() || !dotmuExogenous.empty()){
        if(muExogenous.empty() || dotmuExogenous.empty()){
            throw std::runtime_error("muExogenous and dotmuExogenous must be specified together");
        }
        for(unsigned j = 0; j < ydim; j++) {
            covAllDimensions[j].mu = muExogenous.col(j);
            covAllDimensions[j].dotmu = dotmuExogenous.col(j);
        }
    }
}

void MagiSolver::initTheta() {
    arma::vec sigmaUsed(ydim);
    if(useScalerSigma){
        sigmaUsed.fill(sigmaInit(0));
    }else{
        sigmaUsed = sigmaInit;
    }

    if (!thetaInitExogenous.empty()){
        thetaInit = thetaInitExogenous;
    }else{
        thetaInit = optimizeThetaInit(yFull,
                                      odeModel,
                                      covAllDimensions,
                                      sigmaUsed,
                                      priorTemperature,
                                      xInit,
                                      useBand);
    }
}

void MagiSolver::initMissingComponent() {
    const unsigned int nSGD = 0;  // skip sgd, not useful, and produce unstable result due to delay eval of arma
    double learningRate = 1e-6;
    const arma::uvec & nobsEachDim = arma::sum(indicatorMatWithObs, 0).t();
    const arma::uvec & missingComponentDim = arma::find(nobsEachDim < 3);
    if(missingComponentDim.empty()){
        return;
    }

    const arma::uvec & observedComponentDim = arma::find(nobsEachDim >= 0);
    if(xInitExogenous.empty()) {
        for (auto iPtr = missingComponentDim.begin(); iPtr < missingComponentDim.end(); iPtr++) {
            xInit.col(*iPtr) = arma::mean(xInit.cols(observedComponentDim), 1);
            xInit.submat(idxColElemWithObs[*iPtr], arma::uvec({*iPtr})) =
                    yFull.submat(idxColElemWithObs[*iPtr], arma::uvec({*iPtr}));
        }
    }

    // phi for missing component
    if(phiExogenous.empty()){
        const arma::mat & phiMissingDimensions = optimizePhi(yFull,
                                                             tvecFull,
                                                             odeModel,
                                                             sigmaInit,
                                                             priorTemperature,
                                                             xInit,
                                                             thetaInit,
                                                             phiAllDimensions,
                                                             missingComponentDim);
        if(verbose){
            std::cout << "initMissingComponent: phiMissingDimensions = \n"
                      << phiMissingDimensions << "\n";
        }

        phiAllDimensions.cols(missingComponentDim) = phiMissingDimensions;
    }

    if(verbose) {
        std::cout << "phiAllDimensions = \n" << phiAllDimensions << "\n";
    }

    for(unsigned i = 0; i < missingComponentDim.size(); i++){
        unsigned j = missingComponentDim[i];
        auto mu = covAllDimensions[j].mu;
        auto dotmu = covAllDimensions[j].dotmu;
        covAllDimensions[j] = kernelCov(phiAllDimensions.col(j), distSignedFull.slice(j), 3);
        covAllDimensions[j].addBandCov(bandSize);
        covAllDimensions[j].mu = mu;
        covAllDimensions[j].dotmu = dotmu;
        covAllDimensions[j].tvecCovInput = tvecFull;
    }

    // update theta
    initTheta();

    // x for missing component
    lp llikOld = xthetaphisigmallik( xInit,
                                     thetaInit,
                                     phiAllDimensions,
                                     sigmaInit,
                                     yFull,
                                     tvecFull,
                                     odeModel);
    arma::mat xInitOld = xInit;
    arma::mat thetaInitOld = thetaInit;
    arma::mat phiAllDimensionsOld = phiAllDimensions;

    std::stringstream ss;

    for(unsigned int iSGD = 0; iSGD < nSGD; iSGD++){
        // this std stringstream is essential for the code to run
        std::cout << "initMissingComponent: xInit.row(0) = " << xInit.row(0)
                  << "thetaInit = " << thetaInit.t();
        const lp & llik = xthetaphisigmallik( xInit,
                                              thetaInit,
                                              phiAllDimensions,
                                              sigmaInit,
                                              yFull,
                                              tvecFull,
                                              odeModel);

        if(llik.value - llikOld.value < -std::abs(llikOld.value) * 0.1){
            learningRate *= 0.1;
            xInit = xInitOld;
            thetaInit = thetaInitOld;
            phiAllDimensions = phiAllDimensionsOld;

            if (verbose) {
                std::cout << "initMissingComponent iteration " << iSGD << "; roll back:\n"
                          << "llik.value = " << llik.value
                          << "; learningRate = " << learningRate << "\n";
            }

            continue;
        }

        xInitOld = xInit;
        thetaInitOld = thetaInit;
        phiAllDimensionsOld = phiAllDimensions;
        llikOld = llik;

        if(learningRate < 1e-14){
            break;
        }
        if(iSGD % 100 == 0){
            learningRate *= 10.0;
        }
        thetaInit += learningRate * llik.gradient.subvec(xInit.size(), xInit.size() + thetaInit.size() - 1);
        for (auto iPtr = missingComponentDim.begin(); iPtr < missingComponentDim.end(); iPtr++){
            xInit.col(*iPtr) += learningRate * llik.gradient.subvec(
                    xInit.n_rows * (*iPtr), xInit.n_rows * (*iPtr + 1) - 1);
            phiAllDimensions.col(*iPtr) += learningRate * llik.gradient.subvec(
                    xInit.size() + thetaInit.size() + phiAllDimensions.n_rows * (*iPtr),
                    xInit.size() + thetaInit.size() + phiAllDimensions.n_rows * (*iPtr) + 1);
        }
        thetaInit = arma::max(thetaInit, odeModel.thetaLowerBound + 1e-6);
        thetaInit = arma::min(thetaInit, odeModel.thetaUpperBound - 1e-6);
        phiAllDimensions = arma::max(phiAllDimensions, arma::ones(arma::size(phiAllDimensions)) * 1e-6);

        if (verbose) {
            std::cout << "\ninitMissingComponent iteration " << iSGD
                      << "; xthetaphisigmallik = " << llik.value
                      << "; phi missing dim = \n" << phiAllDimensions.cols(missingComponentDim).t()
                      << "\n";
        }
    }

    if(!skipMissingComponentOptimization){
        try {
            const arma::vec & xthetaphi = optimizeXmissingThetaPhi(yFull,
                                                                   tvecFull,
                                                                   odeModel,
                                                                   sigmaInit,
                                                                   priorTemperature,
                                                                   xInit,
                                                                   thetaInit,
                                                                   phiAllDimensions,
                                                                   missingComponentDim);
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
            }
        } catch (...) {
            std::cout << "Exception occurred in joint optimization for Xmissing,Theta,Phi";
        }
    }

    const lp & llik = xthetaphisigmallik( xInit,
                                          thetaInit,
                                          phiAllDimensions,
                                          sigmaInit,
                                          yFull,
                                          tvecFull,
                                          odeModel);

    std::cout << "\nafter optimization "
              << "; xthetaphisigmallik = " << llik.value
              << "; phi missing dim = \n" << phiAllDimensions.cols(missingComponentDim).t()
              << "\n";
}

void MagiSolver::doHMC(int iEpoch) {
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
    hmcSampler.sampleChian(xthetasigmaInit, stepLow, verbose);
    llikxthetasigmaSamples(arma::span(0, 0), arma::span::all, arma::span(iEpoch, iEpoch)) = hmcSampler.lliklist;
    llikxthetasigmaSamples(arma::span(1, llikxthetasigmaSamples.n_rows - 1), arma::span::all, arma::span(iEpoch, iEpoch)) = hmcSampler.xth;
    stepLow = hmcSampler.stepLow;
}

void MagiSolver::sampleInEpochs() {
    std::string epochMethod = "mean";

    stepLow = arma::vec(llikxthetasigmaSamples.n_rows - 1);
    if (stepSizeFactorHmc.n_elem > 1){
        stepLow = (1.0 / nstepsHmc * stepSizeFactorHmc);
    }else{
        stepLow.fill(1.0 / nstepsHmc * stepSizeFactorHmc(0));
    }

    if(useFixedSigma){
        stepLow.subvec(xInit.size() + thetaInit.size(), stepLow.size() - 1).fill(0);
    }

    for(int iEpoch = 0; iEpoch < nEpoch; iEpoch++){
        doHMC(iEpoch);
        const arma::mat & xthetasigmaSamples = llikxthetasigmaSamples(
                arma::span(1, llikxthetasigmaSamples.n_rows - 1), arma::span::all, arma::span(iEpoch, iEpoch));
        // update mu and dotmu
        arma::mat xPosteriorMean = arma::mean(
                xthetasigmaSamples(
                        arma::span(0, yFull.size() - 1),
                        arma::span(static_cast<int>(niterHmc * burninRatioHmc), niterHmc - 1)),
                1);
        xPosteriorMean.reshape(yFull.n_rows, yFull.n_cols);
        arma::vec thetaPosteriorMean = arma::mean(
                xthetasigmaSamples(
                        arma::span(yFull.size(), yFull.size() + thetaInit.size() - 1),
                        arma::span(static_cast<int>(niterHmc * burninRatioHmc), niterHmc - 1)),
                1);
        arma::vec sigmaPosteriorMean = arma::mean(
                xthetasigmaSamples(
                        arma::span(yFull.size() + thetaInit.size(), xthetasigmaSamples.n_rows - 1),
                        arma::span(static_cast<int>(niterHmc * burninRatioHmc), niterHmc - 1)),
                1);

        // TODO allow median or numerical solver
        for(unsigned long j = 0; j < covAllDimensions.size(); j++){
            covAllDimensions[j].mu = xPosteriorMean.col(j);
        }

        if (epochMethod == "mean"){
            for(unsigned long j = 0; j < covAllDimensions.size(); j++){
                covAllDimensions[j].dotmu = covAllDimensions[j].mphi * xPosteriorMean.col(j);
            }
        }else if(epochMethod == "f_bar_x") {
            arma::mat dotxOde = odeModel.fOde(thetaPosteriorMean, xPosteriorMean, tvecFull);
            for(unsigned long j = 0; j < covAllDimensions.size(); j++) {
                covAllDimensions[j].dotmu = dotxOde.col(j);
            }
        }else if(epochMethod == "bar_f_x"){
            arma::mat dotxOde(yFull.n_rows, yFull.n_cols, arma::fill::zeros);
            for(unsigned it = static_cast<int>(niterHmc * burninRatioHmc); it < niterHmc; it++){
                dotxOde += odeModel.fOde(
                        thetaPosteriorMean,
                        arma::reshape(xthetasigmaSamples(
                                arma::span(0, yFull.size() - 1),
                                arma::span(it, it)
                        ), yFull.n_rows, yFull.n_cols),
                        tvecFull);
            }
            dotxOde /= niterHmc - static_cast<int>(niterHmc * burninRatioHmc);
            for(unsigned long j = 0; j < covAllDimensions.size(); j++) {
                covAllDimensions[j].dotmu = dotxOde.col(j);
            }
        }

        xInit = xPosteriorMean;
        thetaInit = thetaPosteriorMean;
        sigmaInit = sigmaPosteriorMean;
    }
}