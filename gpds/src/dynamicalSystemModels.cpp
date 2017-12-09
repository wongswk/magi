#include "dynamicalSystemModels.h"

using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat fnmodelODE(const arma::vec & theta, const arma::mat & x) {
  const vec & V = x.col(0);
  const vec & R = x.col(1);
  
  const vec & Vdt = theta(2) * (V - pow(V,3) / 3.0 + R);
  const vec & Rdt = -1.0/theta(2) * ( V - theta(0) + theta(1) * R);
  
  return join_horiz(Vdt, Rdt);
}

// [[Rcpp::export]]
arma::cube fnmodelDx(const arma::vec & theta, const arma::mat & x) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols);
  
  const vec & V = x.col(0);
  // const vec & R = x.col(1);
  
  resultDx.slice(0).col(0) = theta(2) * (1 - square(V));
  resultDx.slice(0).col(1).fill( theta(2) );
  
  resultDx.slice(1).col(0).fill(-1.0 / theta(2));
  resultDx.slice(1).col(1).fill( -1.0*theta(1)/theta(2) );
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube fnmodelDtheta(const arma::vec & theta, const arma::mat & x) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & V = x.col(0);
  const vec & R = x.col(1);
  
  resultDtheta.slice(0).col(2) = V - pow(V,3) / 3.0 + R;
  
  resultDtheta.slice(1).col(0).fill( 1.0 / theta(2) );
  resultDtheta.slice(1).col(1) = -R / theta(2);
  resultDtheta.slice(1).col(2) = 1.0/pow(theta(2), 2) * ( V - theta(0) + theta(1) * R);
    
  return resultDtheta;
}

// [[Rcpp::export]]
arma::mat hes1modelODE(const arma::vec & theta, const arma::mat & x) {
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  mat PMHdt(x.n_rows, x.n_cols);
  PMHdt.col(0) = -theta(0)*P%H + theta(1)*M - theta(2)*P;
  PMHdt.col(1) = theta(3)*M + theta(4)/(1+square(P));
  PMHdt.col(2) = -theta(0)*P%H + theta(5)/(1+square(P)) - theta(6)*H;
  
  return PMHdt;
}

// [[Rcpp::export]]
arma::cube hes1modelDx(const arma::vec & theta, const arma::mat & x) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & H = x.col(2); 
  
  resultDx.slice(0).col(0) = -theta(0)*H - theta(2);
  resultDx.slice(0).col(1).fill( theta(1) );
  resultDx.slice(0).col(2) = -theta(0)*P;
  
  resultDx.slice(1).col(0) = -2*theta(4)*P / square(1.0 + square(P));
  resultDx.slice(1).col(1).fill( theta(3) );
  
  resultDx.slice(2).col(0) = -theta(0)*H - 2*theta(5)*P / square(1.0 + square(P));
  resultDx.slice(2).col(2) = -theta(0)*P - theta(6);
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube hes1modelDtheta(const arma::vec & theta, const arma::mat & x) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  resultDtheta.slice(0).col(0) = -P % H;
  resultDtheta.slice(0).col(1) = M;
  resultDtheta.slice(0).col(2) = -P;
  
  resultDtheta.slice(1).col(3) = M;
  resultDtheta.slice(1).col(4) = 1/(1 + square(P));
  
  resultDtheta.slice(2).col(0) = -P % H;
  resultDtheta.slice(2).col(5) = 1/(1 + square(P));
  resultDtheta.slice(2).col(6) = -H;
  
  return resultDtheta;
}
