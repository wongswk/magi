#include "classDefinition.h"


lp xthetaphisigmallik( const arma::mat & xlatent, 
                       const arma::vec & theta, 
                       const arma::mat & phi, 
                       const arma::vec & sigma, 
                       const arma::mat & yobs, 
                       const arma::vec & xtimes,
                       const OdeSystem & fOdeModel);
