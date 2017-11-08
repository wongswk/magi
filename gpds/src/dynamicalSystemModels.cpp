#include "dynamicalSystemModels.h"

using namespace arma;
using namespace std;

mat fnmodelODE(const vec & theta, const mat & x) {
  const vec & V = x.col(0);
  const vec & R = x.col(1);
  
  const vec & Vdt = theta(2) * (V - pow(V,3) / 3.0 + R);
  const vec & Rdt = -1.0/theta(2) * ( V - theta(0) + theta(1) * R);
  
  return join_horiz(Vdt, Rdt);
}

cube fnmodelDx(const vec & theta, const mat & x) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols);
  
  const vec & V = x.col(0);
  // const vec & R = x.col(1);
  
  resultDx.slice(0).col(0) = theta(2) * (1 - square(V));
  resultDx.slice(0).col(1).fill( theta(2) );
  
  resultDx.slice(1).col(0).fill(-1.0 / theta(2));
  resultDx.slice(1).col(0).fill( -1.0*theta(1)/theta(2) );
  
  return resultDx;
}

cube fnmodelDtheta(const vec & theta, const mat & x) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & V = x.col(0);
  const vec & R = x.col(1);
  
  resultDtheta.slice(0).col(2) = V - pow(V,3) / 3.0 + R;
  
  resultDtheta.slice(1).col(0).fill( 1.0 / theta(2) );
  resultDtheta.slice(1).col(1) = -R / theta(2);
  resultDtheta.slice(1).col(2) = 1.0/pow(theta(2), 2) * ( V - theta(0) + theta(1) * R);
    
  return resultDtheta;
}
