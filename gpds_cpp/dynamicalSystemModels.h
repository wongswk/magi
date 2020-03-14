#ifndef DYNAMICALSYSTEMMODELS_H
#define DYNAMICALSYSTEMMODELS_H

#include "classDefinition.h"

arma::mat fnmodelODE(const arma::vec &, const arma::mat &);
arma::cube fnmodelDx(const arma::vec &, const arma::mat &);
arma::cube fnmodelDtheta(const arma::vec &, const arma::mat &);

arma::mat hes1modelODE(const arma::vec &, const arma::mat &);
arma::cube hes1modelDx(const arma::vec &, const arma::mat &);
arma::cube hes1modelDtheta(const arma::vec &, const arma::mat &);

arma::mat hes1logmodelODE(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDx(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDtheta(const arma::vec &, const arma::mat &);

arma::mat hes1logmodelODEfixg(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDxfixg(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDthetafixg(const arma::vec &, const arma::mat &);

arma::mat hes1logmodelODEfixf(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDxfixf(const arma::vec &, const arma::mat &);
arma::cube hes1logmodelDthetafixf(const arma::vec &, const arma::mat &);

arma::mat HIVmodelODE(const arma::vec &, const arma::mat &);
arma::cube HIVmodelDx(const arma::vec &, const arma::mat &);
arma::cube HIVmodelDtheta(const arma::vec &, const arma::mat &);

arma::mat ptransmodelODE(const arma::vec &, const arma::mat &);
arma::cube ptransmodelDx(const arma::vec &, const arma::mat &);
arma::cube ptransmodelDtheta(const arma::vec &, const arma::mat &);

#endif