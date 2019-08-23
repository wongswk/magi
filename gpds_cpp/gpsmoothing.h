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
                            const arma::mat & xInitInput);

arma::vec gpsmooth(const arma::mat & yobsInput,
                   const arma::mat & distInput,
                   std::string kernelInput,
                   bool useFrequencyBasedPrior = false);

arma::cube calcMeanCurve(const arma::vec & xInput,
                         const arma::vec & yInput,
                         const arma::vec & xOutput,
                         const arma::mat & phiCandidates,
                         const arma::vec & sigmaCandidates,
                         const std::string kerneltype = "generalMatern",
                         const bool useDeriv = false);

#endif //GPDS_MULTI_LANG_GPSMOOTHING_H
