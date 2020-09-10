//
// Created by Shihao Yang on 6/5/20.
//

#include "magi_main_py.h"

MagiSolver solveMagiPy(const arma::mat & yFull,
                      const OdeSystem & odeModel,
                      const arma::vec & tvecFull,
                      const arma::vec & sigmaExogenous ,
                      const arma::mat & phiExogenous ,
                      const arma::mat & xInitExogenous ,
                      const arma::vec & thetaInitExogenous ,
                      const arma::mat & muExogenous ,
                      const arma::mat & dotmuExogenous ,
                      const double priorTemperatureLevel ,
                      const double priorTemperatureDeriv ,
                      const double priorTemperatureObs ,
                      std::string kernel ,
                      const int nstepsHmc ,
                      const double burninRatioHmc ,
                      const unsigned int niterHmc ,
                      const double stepSizeFactorHmc ,
                      const int nEpoch ,
                      const int bandSize ,
                      bool useFrequencyBasedPrior ,
                      bool useBand ,
                      bool useMean ,
                      bool useScalerSigma ,
                      bool useFixedSigma ,
                      bool verbose) {

    MagiSolver solver(yFull,
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
//    return solver.llikxthetasigmaSamples;
    return solver;
}
