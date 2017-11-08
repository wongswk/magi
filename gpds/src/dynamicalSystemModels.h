#ifndef DYNAMICALSYSTEMMODELS_H
#define DYNAMICALSYSTEMMODELS_H

#include "classDefinition.h"

arma::mat fnmodelODE(const arma::vec &, const arma::mat &);
arma::cube fnmodelDx(const arma::vec &, const arma::mat &);
arma::cube fnmodelDtheta(const arma::vec &, const arma::mat &);

#endif