#include "classDefinition.h"
#include "GpdsSolver.h"


// [[Rcpp::export]]
arma::cube solveGpds(const arma::mat & yFull,
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
                    const double stepSizeFactorHmc = 1,
                    const int nEpoch = 10,
                    const int bandSize = 20,
                    bool useFrequencyBasedPrior = false,
                    bool useBand = true,
                    bool useMean = true,
                    bool useScalerSigma = false,
                    bool useFixedSigma = false,
                    bool verbose = false) {

    GpdsSolver solver(yFull,
                      odeModel,
                      tvecFull,
                      sigmaExogenous,
                      phiExogenous,
                      xInitExogenous,
                      thetaInitExogenous,
                      muExogenous,
                      dotmuExogenous,
                      priorTemperatureLevel,
                      priorTemperatureDeriv,
                      priorTemperatureObs,
                      std::move(kernel),
                      nstepsHmc,
                      burninRatioHmc,
                      niterHmc,
                      stepSizeFactorHmc,
                      nEpoch,
                      bandSize,
                      useFrequencyBasedPrior,
                      useBand,
                      useMean,
                      useScalerSigma,
                      useFixedSigma,
                      verbose);
    solver.setupPhiSigma();
    if(verbose){
        std::cout << "phi = \n" << solver.phiAllDimensions << "\n";
    }
    solver.initXmudotmu();
    solver.initTheta();
    solver.initMissingComponent();
    solver.sampleInEpochs();
    return solver.llikxthetasigmaSamples;
}
