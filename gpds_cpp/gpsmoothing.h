//
// Created by Shihao Yang on 2019-08-22.
//

#ifndef GPDS_MULTI_LANG_GPSMOOTHING_H
#define GPDS_MULTI_LANG_GPSMOOTHING_H

#include "classDefinition.h"
arma::vec optimizeThetaInit(const arma::mat & yobsInput,
                            const OdeSystem & fOdeModelInput,
                            const std::vector<gpcov> & covAllDimensionsInput,
                            const arma::vec & sigmaAllDimensionsInput,
                            const arma::vec & priorTemperatureInput,
                            const arma::mat & xInitInput,
                            const bool useBandInput);

arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   const double sigmaExogenScalar = -1,
                   bool useFrequencyBasedPrior = false);

arma::cube calcMeanCurve(const arma::vec & xInput,
                         const arma::vec & yInput,
                         const arma::vec & xOutput,
                         const arma::mat & phiCandidates,
                         const arma::vec & sigmaCandidates,
                         const std::string kerneltype = "generalMatern",
                         const bool useDeriv = false);

arma::mat optimizePhi(const arma::mat & yobsInput,
                      const arma::vec & tvecInput,
                      const OdeSystem & fOdeModelInput,
                      const arma::vec & sigmaAllDimensionsInput,
                      const arma::vec & priorTemperatureInput,
                      const arma::mat & xInitInput,
                      const arma::vec & thetaInitInput,
                      const arma::mat & phiInitInput,
                      const arma::uvec & missingComponentDim);

arma::mat optimizeXmissingThetaPhi(const arma::mat & yobsInput,
                                   const arma::vec & tvecInput,
                                   const OdeSystem & fOdeModelInput,
                                   const arma::vec & sigmaAllDimensionsInput,
                                   const arma::vec & priorTemperatureInput,
                                   const arma::mat & xInitInput,
                                   const arma::vec & thetaInitInput,
                                   const arma::mat & phiInitInput,
                                   const arma::uvec & missingComponentDim);

#endif //GPDS_MULTI_LANG_GPSMOOTHING_H
