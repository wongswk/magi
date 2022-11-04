#ifndef MAGI_MULTI_LANG_MAGISOLVER_H
#define MAGI_MULTI_LANG_MAGISOLVER_H

#include "classDefinition.h"

class MagiSolver {
public:
    // input data
    const arma::mat & yFull;
    const OdeSystem & odeModel;
    const arma::vec & tvecFull;

    // optional exogenous value
    const arma::vec & sigmaExogenous;
    const arma::mat & phiExogenous;
    const arma::mat & xInitExogenous;
    const arma::vec & thetaInitExogenous;
    const arma::mat & muExogenous;
    const arma::mat & dotmuExogenous;

    // configuration
    const arma::vec priorTemperature;
    std::string kernel;
    const int nstepsHmc;
    const double burninRatioHmc;
    const unsigned int niterHmc;
    arma::vec stepSizeFactorHmc;
    const int nEpoch;
    const int bandSize;
    bool useFrequencyBasedPrior;
    bool useBand;
    bool useMean;
    bool useScalerSigma;
    bool useFixedSigma;
    bool verbose;
    bool skipMissingComponentOptimization;

    // intermediate object storage
    const unsigned int ydim;
    const unsigned int sigmaSize;
    std::vector<gpcov> covAllDimensions;
    std::string loglikflag;
    arma::cube distSignedFull;
    std::function<gpcov(arma::vec, arma::mat, int)> kernelCov;

    arma::mat yObs;
    arma::mat distSignedObs;
    arma::uvec indicatorRowWithObs;
    arma::uvec idxRowWithObs;
    arma::umat indicatorMatWithObs;
    std::vector<arma::uvec> idxColElemWithObs;


    arma::mat phiAllDimensions;
    arma::vec sigmaInit;
    arma::mat xInit;
    arma::mat thetaInit;

    arma::vec stepLow;

    // output
    arma::cube llikxthetasigmaSamples;

    MagiSolver(const arma::mat & yFull,
               const OdeSystem & odeModel,
               const arma::vec & tvecFull,
               const arma::vec & sigmaExogenous = arma::vec(),
               const arma::mat & phiExogenous = arma::mat(),
               const arma::mat & xInitExogenous = arma::mat(),
               const arma::vec & thetaInitExogenous = arma::vec(),
               const arma::mat & muExogenous = arma::mat(),
               const arma::mat & dotmuExogenous = arma::mat(),
               const double priorTemperatureLevel = 1,
               const double priorTemperatureDeriv = 1,
               const double priorTemperatureObs = 1,
               std::string kernel = "generalMatern",
               const int nstepsHmc = 500,
               const double burninRatioHmc = 0.5,
               const unsigned int niterHmc = 10000,
               const arma::vec stepSizeFactorHmc = arma::vec(),
               const arma::cube distSignedFullInput = arma::cube(),
               const int nEpoch = 1,
               const int bandSize = 20,
               bool useFrequencyBasedPrior = false,
               bool useBand = true,
               bool useMean = true,
               bool useScalerSigma = false,
               bool useFixedSigma = false,
               bool skipMissingComponentOptimization=false,
               bool verbose = false);

    void setupPhiSigma();
    void initXmudotmu();
    void initTheta();
    void initMissingComponent();
    void doHMC(int iEpoch);
    void sampleInEpochs();
};


#endif //MAGI_MULTI_LANG_MAGISOLVER_H
