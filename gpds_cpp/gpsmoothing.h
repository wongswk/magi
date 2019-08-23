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

#endif //GPDS_MULTI_LANG_GPSMOOTHING_H
