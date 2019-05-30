#include "classDefinition.h"

using namespace arma;
using namespace std;

bool OdeSystem::checkBound(const arma::mat & xlatent, const arma::vec & theta, lp* retPtr) const{
  uvec x2Large;
  if( xUpperBound.size() > 0 && any( xUpperBound < datum::inf)){
    x2Large = find(xlatent > xUpperBound);  // FIXME bug here for comparison
  }
  uvec x2Small;
  if( xLowerBound.size() > 0 && any( xLowerBound > -datum::inf)){
    x2Small = find(xlatent < xLowerBound); // FIXME bug here for comparison
  }
  uvec theta2Large = find(theta > thetaUpperBound);
  uvec theta2Small = find(theta < thetaLowerBound);
  
  if(x2Large.size() || x2Small.size() || theta2Large.size() || theta2Small.size()){
    retPtr->value = -1e+9;
    retPtr->gradient = arma::zeros(xlatent.size()+theta.size());
    theta2Large = theta2Large + xlatent.size();
    theta2Small = theta2Small + xlatent.size();
    retPtr->gradient(x2Large).fill(-1e+9);
    retPtr->gradient(x2Small).fill( 1e+9);
    retPtr->gradient(theta2Large).fill(-1e+9);
    retPtr->gradient(theta2Small).fill( 1e+9);
    return true;
  }
  return false;
}