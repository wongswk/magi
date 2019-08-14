//
// Created by Shihao Yang on 8/14/19.
//

#ifndef DYNAMIC_SYSTEMS_PHI1LOGLIKELIHOOD_H
#define DYNAMIC_SYSTEMS_PHI1LOGLIKELIHOOD_H

#include "classDefinition.h"

lp xthetaphi1sigmallik( const arma::mat & xlatent,
                        const arma::vec & theta,
                        const arma::vec & phi1,
                        const arma::vec & sigmaInput,
                        const arma::mat & yobs,
                        const std::vector<gpcov> & CovAllDimensions,
                        const OdeSystem & fOdeModel,
                        const arma::vec & priorTemperatureInput = arma::ones(1),
                        const bool useBand = false,
                        const bool useMean = false);

#endif //DYNAMIC_SYSTEMS_PHI1LOGLIKELIHOOD_H
