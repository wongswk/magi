#ifndef DYNAMIC_SYSTEMS_RCPPTRANSLATION_H
#define DYNAMIC_SYSTEMS_RCPPTRANSLATION_H

#include "RcppArmadillo.h"
#include "classDefinition.h"

namespace Rcpp
{
    // gpcov
    template <>
    gpcov as(SEXP x);

    template <>
    SEXP wrap(const gpcov& object);

}


#endif //DYNAMIC_SYSTEMS_RCPPTRANSLATION_H
