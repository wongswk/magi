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

    template <>
    std::vector<gpcov> as(SEXP x);

    // lp
    template <>
    lp as(SEXP x);

    // OdeSystem
    template <>
    OdeSystem as(SEXP x);

}


#endif //DYNAMIC_SYSTEMS_RCPPTRANSLATION_H
