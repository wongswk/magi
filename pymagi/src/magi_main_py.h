//
// Created by Shihao Yang on 6/5/20.
//

#ifndef MAGI_MULTI_LANG_MAGI_MAIN_PY_H
#include <classDefinition.h>
#include <MagiSolver.h>
#define MAGI_MULTI_LANG_MAGI_MAIN_PY_H

MagiSolver solveMagiPy(const arma::mat & yFull,
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
                       const arma::vec & stepSizeFactorHmc = arma::vec(),
                       const int nEpoch = 10,
                       const int bandSize = 20,
                       bool useFrequencyBasedPrior = false,
                       bool useBand = true,
                       bool useMean = true,
                       bool useScalerSigma = false,
                       bool useFixedSigma = false,
                       bool skipMissingComponentOptimization = false,
                       bool positiveSystem = false,
                       bool verbose = false);

#endif //MAGI_MULTI_LANG_MAGI_MAIN_PY_H
