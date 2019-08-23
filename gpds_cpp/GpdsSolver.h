#ifndef GPDS_MULTI_LANG_GPDSSOLVER_H
#define GPDS_MULTI_LANG_GPDSSOLVER_H

#include "classDefinition.h"

class GpdsSolver {
public:
    // input data
    const arma::mat & yFull;
    const OdeSystem & odeModel;
    const arma::vec & tvecFull;
    const arma::vec & sigmaExogenous;

    // configuration
    const arma::vec & priorTemperature;
    std::string kernel;
    const int nstepsHmc;
    const double burninRatioHmc;
    const unsigned int niterHmc;
    const double stepSizeFactorHmc;
    const int nEpoch;
    const int bandSize;
    bool useFrequencyBasedPrior;
    bool useBand;
    bool useMean;
    bool useScalerSigma;
    bool useFixedSigma;
    bool verbose;

    // intermediate object storage
    const unsigned int ydim;
    const unsigned int sigmaSize;
    std::vector<gpcov> covAllDimensions;
    std::string loglikflag;
    arma::mat distSignedFull;
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
    arma::vec llikSamples;
    arma::mat xthetasigmaSamples;

    GpdsSolver(const arma::mat & yFull,
               const OdeSystem & odeModel,
               const arma::vec & tvecFull,
               const arma::vec & sigmaExogenous = arma::vec(),
               const double priorTemperatureLevel = 1,
               const double priorTemperatureDeriv = 1,
               std::string kernel = "generalMatern",
               const int nstepsHmc = 500,
               const double burninRatioHmc = 0.5,
               const unsigned int niterHmc = 10000,
               const double stepSizeFactorHmc = 1,
               const int nEpoch = 10,
               const int bandSize = 20,
               bool useFrequencyBasedPrior = false,
               bool useBand = true,
               bool useMean = true,
               bool useScalerSigma = false,
               bool useFixedSigma = false,
               bool verbose = false);

    void setupPhiSigma();
    void initXmudotmu();
    void initTheta();
    void doHMC();
};


#endif //GPDS_MULTI_LANG_GPDSSOLVER_H
