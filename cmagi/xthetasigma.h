//
// Created by Shihao Yang on 8/14/19.
//

#ifndef DYNAMIC_SYSTEMS_XTHETASIGMA_H

#include "classDefinition.h"

//' log likelihood for latent states and ODE theta conditional on phi sigma
//'
//' @param phisig      the parameter phi and sigma
//' @param yobs        observed data
lp xthetasigmallik( const arma::mat & xlatent,
                    const arma::vec & theta,
                    const arma::vec & sigmaInput,
                    const arma::mat & yobs,
                    const std::vector<gpcov> & CovAllDimensions,
                    const OdeSystem & fOdeModel,
                    const arma::vec & priorTemperatureInput = arma::ones(1),
                    const bool useBand = false,
                    const bool useMean = false);

#define DYNAMIC_SYSTEMS_XTHETASIGMA_H

#endif //DYNAMIC_SYSTEMS_XTHETASIGMA_H
