#ifndef DYNAMIC_SYSTEMS_RCPPTRANSLATION_H
#define DYNAMIC_SYSTEMS_RCPPTRANSLATION_H

#include "RcppArmadillo.h"
#include "classDefinition.h"
#include "Sampler.h"

namespace Rcpp
{
    // gpcov
    template <>
    gpcov as(SEXP x);

    template <>
    SEXP wrap(const gpcov& object);

    // lp
    template <>
    lp as(SEXP x);

    // OdeSystem
    template <>
    OdeSystem as(SEXP x);

    // Sampler
    template <>
    SEXP wrap(const Sampler& object);

}


#endif //DYNAMIC_SYSTEMS_RCPPTRANSLATION_H
