#include "dynamicalSystemModels.h"

using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat fnmodelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & V = x.col(0);
  const vec & R = x.col(1);
  
  const vec & Vdt = theta(2) * (V - pow(V,3) / 3.0 + R);
  const vec & Rdt = -1.0/theta(2) * ( V - theta(0) + theta(1) * R);
  
  return join_horiz(Vdt, Rdt);
}

// [[Rcpp::export]]
arma::cube fnmodelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
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
arma::cube fnmodelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
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
arma::mat hes1modelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  mat PMHdt(x.n_rows, x.n_cols);
  PMHdt.col(0) = -theta(0)*P%H + theta(1)*M - theta(2)*P;
  PMHdt.col(1) = -theta(3)*M + theta(4)/(1+square(P));
  PMHdt.col(2) = -theta(0)*P%H + theta(5)/(1+square(P)) - theta(6)*H;
  
  return PMHdt;
}

// [[Rcpp::export]]
arma::cube hes1modelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & H = x.col(2); 
  
  resultDx.slice(0).col(0) = -theta(0)*H - theta(2);
  resultDx.slice(0).col(1).fill( theta(1) );
  resultDx.slice(0).col(2) = -theta(0)*P;
  
  resultDx.slice(1).col(0) = -2*theta(4)*P / square(1.0 + square(P));
  resultDx.slice(1).col(1).fill( -theta(3) );
  
  resultDx.slice(2).col(0) = -theta(0)*H - 2*theta(5)*P / square(1.0 + square(P));
  resultDx.slice(2).col(2) = -theta(0)*P - theta(6);
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube hes1modelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  resultDtheta.slice(0).col(0) = -P % H;
  resultDtheta.slice(0).col(1) = M;
  resultDtheta.slice(0).col(2) = -P;
  
  resultDtheta.slice(1).col(3) = -M;
  resultDtheta.slice(1).col(4) = 1/(1 + square(P));
  
  resultDtheta.slice(2).col(0) = -P % H;
  resultDtheta.slice(2).col(5) = 1/(1 + square(P));
  resultDtheta.slice(2).col(6) = -H;
  
  return resultDtheta;
}

// [[Rcpp::export]]
arma::mat hes1logmodelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & P = arma::exp(x.col(0));
  const vec & M = arma::exp(x.col(1));
  const vec & H = arma::exp(x.col(2)); 
  
  mat PMHdt(x.n_rows, x.n_cols);
  PMHdt.col(0) = -theta(0)*H + theta(1)*M/P - theta(2);
  PMHdt.col(1) = -theta(3) + theta(4)/(1+square(P))/M;
  PMHdt.col(2) = -theta(0)*P + theta(5)/(1+square(P))/H - theta(6);
  
  return PMHdt;
}

// [[Rcpp::export]]
arma::cube hes1logmodelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1); 
  const vec & H = x.col(2); 
  
  const vec & expMminusP = exp(M-P);
  const vec & dP = -pow(1+exp(2*P), -2)%exp(2*P)*2;
  
  resultDx.slice(0).col(0) = -theta(1)*expMminusP;
  resultDx.slice(0).col(1) = theta(1)*expMminusP;
  resultDx.slice(0).col(2) = -theta(0)*exp(H);
  
  resultDx.slice(1).col(0) = theta(4)*exp(-M)%dP;
  resultDx.slice(1).col(1) = -theta(4)*exp(-M)/(1+exp(2*P));
  
  resultDx.slice(2).col(0) = -theta(0)*exp(P) + theta(5)*exp(-H)%dP;
  resultDx.slice(2).col(2) = -theta(5)*exp(-H)/(1+exp(2*P));
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube hes1logmodelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  resultDtheta.slice(0).col(0) = -exp(H);
  resultDtheta.slice(0).col(1) = exp(M-P);
  resultDtheta.slice(0).col(2).fill(-1);
  
  resultDtheta.slice(1).col(3).fill(-1);
  resultDtheta.slice(1).col(4) = exp(-M)/(1+exp(2*P));
  
  resultDtheta.slice(2).col(0) = -exp(P);
  resultDtheta.slice(2).col(5) = exp(-H)/(1+exp(2*P));
  resultDtheta.slice(2).col(6).fill(-1);
  
  return resultDtheta;
}

// [[Rcpp::export]]
arma::mat hes1logmodelODEfixg(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & P = arma::exp(x.col(0));
  const vec & M = arma::exp(x.col(1));
  const vec & H = arma::exp(x.col(2)); 
  
  mat PMHdt(x.n_rows, x.n_cols);
  PMHdt.col(0) = -theta(0)*H + theta(1)*M/P - theta(2);
  PMHdt.col(1) = -theta(3) + theta(4)/(1+square(P))/M;
  PMHdt.col(2) = -theta(0)*P + theta(5)/(1+square(P))/H - 0.3;
  
  return PMHdt;
}

// [[Rcpp::export]]
arma::cube hes1logmodelDxfixg(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1); 
  const vec & H = x.col(2); 
  
  const vec & expMminusP = exp(M-P);
  const vec & dP = -pow(1+exp(2*P), -2)%exp(2*P)*2;
  
  resultDx.slice(0).col(0) = -theta(1)*expMminusP;
  resultDx.slice(0).col(1) = theta(1)*expMminusP;
  resultDx.slice(0).col(2) = -theta(0)*exp(H);
  
  resultDx.slice(1).col(0) = theta(4)*exp(-M)%dP;
  resultDx.slice(1).col(1) = -theta(4)*exp(-M)/(1+exp(2*P));
  
  resultDx.slice(2).col(0) = -theta(0)*exp(P) + theta(5)*exp(-H)%dP;
  resultDx.slice(2).col(2) = -theta(5)*exp(-H)/(1+exp(2*P));
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube hes1logmodelDthetafixg(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & P = x.col(0);
  const vec & M = x.col(1);
  const vec & H = x.col(2); 
  
  resultDtheta.slice(0).col(0) = -exp(H);
  resultDtheta.slice(0).col(1) = exp(M-P);
  resultDtheta.slice(0).col(2).fill(-1);
  
  resultDtheta.slice(1).col(3).fill(-1);
  resultDtheta.slice(1).col(4) = exp(-M)/(1+exp(2*P));
  
  resultDtheta.slice(2).col(0) = -exp(P);
  resultDtheta.slice(2).col(5) = exp(-H)/(1+exp(2*P));
  //resultDtheta.slice(2).col(6).fill(-1);
  
  return resultDtheta;
}


// [[Rcpp::export]]
arma::mat hes1logmodelODEfixf(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & P = arma::exp(x.col(0));
    const vec & M = arma::exp(x.col(1));
    const vec & H = arma::exp(x.col(2));

    mat PMHdt(x.n_rows, x.n_cols);
    PMHdt.col(0) = -theta(0)*H + theta(1)*M/P - theta(2);
    PMHdt.col(1) = -theta(3) + theta(4)/(1+square(P))/M;
    PMHdt.col(2) = -theta(0)*P + 20.0/(1+square(P))/H - theta(5);

    return PMHdt;
}

// [[Rcpp::export]]
arma::cube hes1logmodelDxfixf(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & P = x.col(0);
    const vec & M = x.col(1);
    const vec & H = x.col(2);

    const vec & expMminusP = exp(M-P);
    const vec & dP = -pow(1+exp(2*P), -2)%exp(2*P)*2;

    resultDx.slice(0).col(0) = -theta(1)*expMminusP;
    resultDx.slice(0).col(1) = theta(1)*expMminusP;
    resultDx.slice(0).col(2) = -theta(0)*exp(H);

    resultDx.slice(1).col(0) = theta(4)*exp(-M)%dP;
    resultDx.slice(1).col(1) = -theta(4)*exp(-M)/(1+exp(2*P));

    resultDx.slice(2).col(0) = -theta(0)*exp(P) + 20.0*exp(-H)%dP;
    resultDx.slice(2).col(2) = -20.0*exp(-H)/(1+exp(2*P));

    return resultDx;
}


// [[Rcpp::export]]
arma::cube hes1logmodelDthetafixf(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & P = x.col(0);
    const vec & M = x.col(1);
    const vec & H = x.col(2);

    resultDtheta.slice(0).col(0) = -exp(H);
    resultDtheta.slice(0).col(1) = exp(M-P);
    resultDtheta.slice(0).col(2).fill(-1);

    resultDtheta.slice(1).col(3).fill(-1);
    resultDtheta.slice(1).col(4) = exp(-M)/(1+exp(2*P));

    resultDtheta.slice(2).col(0) = -exp(P);
    resultDtheta.slice(2).col(5).fill(-1);

    return resultDtheta;
}


// [[Rcpp::export]]
arma::mat HIVmodelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & T = exp(x.col(0));
  const vec & Tm = exp(x.col(1));
  const vec & Tw = exp(x.col(2));
  const vec & Tmw = exp(x.col(3));
  // const vec & T = x.col(0);
  // const vec & Tm = x.col(1);
  // const vec & Tw = x.col(2); 
  // const vec & Tmw = x.col(3);
  
    
  mat HIVdt(x.n_rows, x.n_cols);

  HIVdt.col(0) = (theta(0) - 1e-6*theta(1)*Tm - 1e-6*theta(2)*Tw - 1e-6*theta(3)*Tmw);
  HIVdt.col(1) = (theta(6) + 1e-6*theta(1)*T - 1e-6*theta(4)*Tw) + 1e-6*0.25*theta(3)*Tmw%T / Tm;
  HIVdt.col(2) = (theta(7) + 1e-6*theta(2)*T - 1e-6*theta(5)*Tm) + 1e-6*0.25*theta(3)*Tmw%T / Tw;
  HIVdt.col(3) = theta(8) + 0.5*1e-6*theta(3)*T + (1e-6*theta(4)+1e-6*theta(5))*Tw%Tm / Tmw;
  
  return HIVdt;
}

// [[Rcpp::export]]
arma::cube HIVmodelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & T = exp(x.col(0));
  const vec & Tm = exp(x.col(1));
  const vec & Tw = exp(x.col(2));
  const vec & Tmw = exp(x.col(3));

  resultDx.slice(0).col(0).fill(0);
  resultDx.slice(0).col(1) = -1e-6*theta(1)*Tm;
  resultDx.slice(0).col(2) = -1e-6*theta(2)*Tw;
  resultDx.slice(0).col(3) = -1e-6*theta(3)*Tmw;
  
  resultDx.slice(1).col(0) = 1e-6*theta(1)*T + 1e-6*0.25*theta(3)*Tmw%T/Tm;
  resultDx.slice(1).col(1) = -1e-6*0.25*theta(3)*Tmw%T / Tm;
  resultDx.slice(1).col(2) = -1e-6*theta(4)*Tw;
  resultDx.slice(1).col(3) = 0.25*1e-6*theta(3)*Tmw%T/Tm;
  
  resultDx.slice(2).col(0) = 1e-6*theta(2)*T + 0.25*1e-6*theta(3)*Tmw%T/Tw;
  resultDx.slice(2).col(1) = -1e-6*theta(5)*Tm;
  resultDx.slice(2).col(2) = -1e-6*0.25*theta(3)*Tmw%T / Tw;
  resultDx.slice(2).col(3) = 1e-6*0.25*theta(3)*Tmw%T/Tw;
  
  resultDx.slice(3).col(0) = 1e-6*0.5*theta(3)*T;
  resultDx.slice(3).col(1) = (1e-6*theta(4)+1e-6*theta(5))*Tw%Tm/Tmw;
  resultDx.slice(3).col(2) = (1e-6*theta(4)+1e-6*theta(5))*Tm%Tw/Tmw;
  resultDx.slice(3).col(3) = -(1e-6*theta(4)+1e-6*theta(5))*Tw%Tm/Tmw;  
      
  return resultDx;
}

// [[Rcpp::export]]
arma::cube HIVmodelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & T = exp(x.col(0));
  const vec & Tm = exp(x.col(1));
  const vec & Tw = exp(x.col(2));
  const vec & Tmw = exp(x.col(3));
  

  resultDtheta.slice(0).col(0).fill(1.0);
  resultDtheta.slice(0).col(1) = -1e-6*Tm ;
  resultDtheta.slice(0).col(2) = -1e-6*Tw ;
  resultDtheta.slice(0).col(3) = -1e-6*Tmw;
  
  resultDtheta.slice(1).col(1) = 1e-6*T;
  resultDtheta.slice(1).col(3) = 1e-6*0.25*Tmw%T / Tm;
  resultDtheta.slice(1).col(4) = -1e-6*Tw;
  resultDtheta.slice(1).col(6).fill(1.0);

  resultDtheta.slice(2).col(2) = 1e-6*T;
  resultDtheta.slice(2).col(3) = 1e-6*0.25*Tmw%T / Tw;
  resultDtheta.slice(2).col(5) = -1e-6*Tm;
  resultDtheta.slice(2).col(7).fill(1.0);
  
  resultDtheta.slice(3).col(3) = 1e-6*0.5 * T;
  resultDtheta.slice(3).col(4) = 1e-6*Tw % Tm / Tmw;
  resultDtheta.slice(3).col(5) = 1e-6*Tw % Tm / Tmw;
  resultDtheta.slice(3).col(8).fill(1.0);
        
  return resultDtheta;
}

// [[Rcpp::export]]
arma::mat ptransmodelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  const vec & S = x.col(0);
  const vec & dS = x.col(1);
  const vec & R = x.col(2);
  const vec & RS = x.col(3);
  const vec & RPP = x.col(4);
  
  mat resultdt(x.n_rows, x.n_cols);

  resultdt.col(0) = -theta(0)*S - theta(1) * S % R + theta(2) * RS;
  resultdt.col(1) = theta(0)*S;
  resultdt.col(2) = -theta(1)*S%R + theta(2)*RS + theta(4) * RPP / (theta(5)+RPP);
  resultdt.col(3) = theta(1)*S%R - theta(2)* RS - theta(3)*RS;
  resultdt.col(4) = theta(3)*RS - theta(4) * RPP / (theta(5)+RPP);
  
  return resultdt;
}

// [[Rcpp::export]]
arma::cube ptransmodelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);
  
  const vec & S = x.col(0);
  const vec & dS = x.col(1);
  const vec & R = x.col(2);
  const vec & RS = x.col(3);
  const vec & RPP = x.col(4);
  
  resultDx.slice(0).col(0) = -theta(0) - theta(1) * R;
  resultDx.slice(0).col(2) = -theta(1) * S;
  resultDx.slice(0).col(3).fill(theta(2));
  
  resultDx.slice(1).col(0).fill(theta(0));
    
  resultDx.slice(2).col(0) = -theta(1)*R;
  resultDx.slice(2).col(2) = -theta(1)*S;
  resultDx.slice(2).col(3).fill(theta(2));
  resultDx.slice(2).col(4) =  theta(4) * theta(5) /  square(theta(5) + RPP);
  
  resultDx.slice(3).col(0) = theta(1)*R;
  resultDx.slice(3).col(2) = theta(1)*S;
  resultDx.slice(3).col(3).fill(-theta(2) - theta(3));
  
  resultDx.slice(4).col(3).fill(theta(3));
  resultDx.slice(4).col(4) = -theta(4) * theta(5) /  square(theta(5) + RPP);
  
  return resultDx;
}

// [[Rcpp::export]]
arma::cube ptransmodelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
  cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);
  
  const vec & S = x.col(0);
  const vec & dS = x.col(1);
  const vec & R = x.col(2);
  const vec & RS = x.col(3);
  const vec & RPP = x.col(4);  
  
  resultDtheta.slice(0).col(0) = -S;
  resultDtheta.slice(0).col(1) = -S%R;
  resultDtheta.slice(0).col(2) = RS;
  
  resultDtheta.slice(1).col(0) = S;
  
  resultDtheta.slice(2).col(1) = -S%R;
  resultDtheta.slice(2).col(2) = RS;
  resultDtheta.slice(2).col(4) = RPP / (theta(5)+RPP);
  resultDtheta.slice(2).col(5) = -theta(4) * RPP / square(theta(5)+RPP);
  
  resultDtheta.slice(3).col(1) = S%R;
  resultDtheta.slice(3).col(2) = -RS;
  resultDtheta.slice(3).col(3) = -RS;
  
  resultDtheta.slice(4).col(3) = RS;
  resultDtheta.slice(4).col(4) = - RPP / (theta(5)+RPP);
  resultDtheta.slice(4).col(5) = theta(4) * RPP / square(theta(5)+RPP);;
  
  return resultDtheta;
}


// [[Rcpp::export]]
arma::mat MichaelisMentenModelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -theta[0] * e % s + (theta[1]+theta[2]) * es;
    resultdt.col(1) = -theta[0] * e % s + (theta[1]) * es;
    resultdt.col(2) = theta[0] * e % s - (theta[1]+theta[2]) * es;
    resultdt.col(3) = theta[2] * es;

    return resultdt;
}


// [[Rcpp::export]]
arma::cube MichaelisMentenModelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);


    resultDx.slice(0).col(0) = -theta[0] * s;
    resultDx.slice(0).col(1) = -theta[0] * e;
    resultDx.slice(0).col(2).fill(theta[1] + theta[2]);

    resultDx.slice(1).col(0) = -theta[0] * s;
    resultDx.slice(1).col(1) = -theta[0] * e;
    resultDx.slice(1).col(2).fill(theta[1]);

    resultDx.slice(2).col(0) = theta[0] * s;
    resultDx.slice(2).col(1) = theta[0] * e;
    resultDx.slice(2).col(2).fill(-theta[1] - theta[2]);

    resultDx.slice(3).col(2).fill(theta[2]);

    return resultDx;
}


// [[Rcpp::export]]
arma::cube MichaelisMentenModelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);


    resultDtheta.slice(0).col(0) = -e % s;
    resultDtheta.slice(0).col(1) = es;
    resultDtheta.slice(0).col(2) = es;

    resultDtheta.slice(1).col(0) = -e % s;
    resultDtheta.slice(1).col(1) = es;

    resultDtheta.slice(2).col(0) = e % s;
    resultDtheta.slice(2).col(1) = -es;
    resultDtheta.slice(2).col(2) = -es;

    resultDtheta.slice(3).col(2) = es;

    return resultDtheta;
}



// [[Rcpp::export]]
arma::mat MichaelisMentenLogModelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & logE = x.col(0);
    const vec & logS = x.col(1);
    const vec & logES = x.col(2);
    const vec & logP = x.col(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -theta[0] * arma::exp(logS) + (theta[1]+theta[2]) * arma::exp(logES-logE);
    resultdt.col(1) = -theta[0] * arma::exp(logE) + (theta[1]) * arma::exp(logES-logS);
    resultdt.col(2) = theta[0] * arma::exp(logE+logS-logES) - (theta[1]+theta[2]);
    resultdt.col(3) = theta[2] * arma::exp(logES-logP);

    return resultdt;
}


// [[Rcpp::export]]
arma::cube MichaelisMentenLogModelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & logE = x.col(0);
    const vec & logS = x.col(1);
    const vec & logES = x.col(2);
    const vec & logP = x.col(3);


    resultDx.slice(0).col(0) = -(theta[1] + theta[2]) * arma::exp(logES-logE);
    resultDx.slice(0).col(1) = -theta[0] * arma::exp(logS);
    resultDx.slice(0).col(2) = (theta[1] + theta[2]) * arma::exp(logES-logE);

    resultDx.slice(1).col(0) = -theta[0] * arma::exp(logE);
    resultDx.slice(1).col(1) = -theta[1] * arma::exp(logES-logS);
    resultDx.slice(1).col(2) = theta[1] * arma::exp(logES-logS);

    resultDx.slice(2).col(0) = theta[0] * arma::exp(logE+logS-logES);
    resultDx.slice(2).col(1) = theta[0] * arma::exp(logE+logS-logES);
    resultDx.slice(2).col(2) = -theta[0] * arma::exp(logE+logS-logES);

    resultDx.slice(3).col(2) = theta[2] * arma::exp(logES-logP);
    resultDx.slice(3).col(3) = -theta[2] * arma::exp(logES-logP);

    return resultDx;
}


// [[Rcpp::export]]
arma::cube MichaelisMentenLogModelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & logE = x.col(0);
    const vec & logS = x.col(1);
    const vec & logES = x.col(2);
    const vec & logP = x.col(3);


    resultDtheta.slice(0).col(0) = -arma::exp(logS);
    resultDtheta.slice(0).col(1) = arma::exp(logES-logE);
    resultDtheta.slice(0).col(2) = arma::exp(logES-logE);

    resultDtheta.slice(1).col(0) = -arma::exp(logE);
    resultDtheta.slice(1).col(1) = arma::exp(logES-logS);

    resultDtheta.slice(2).col(0) = arma::exp(logE+logS-logES);
    resultDtheta.slice(2).col(1).fill(-1);
    resultDtheta.slice(2).col(2).fill(-1);

    resultDtheta.slice(3).col(2) = arma::exp(logES-logP);

    return resultDtheta;
}

// [[Rcpp::export]]
arma::mat MichaelisMentenModelVaODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -theta[0] * e % s + (theta[1]+theta[2]) * es;
    resultdt.col(1) = -2*theta[0] * e % s + (2*theta[1]) * es;
    resultdt.col(2) = theta[0] * e % s - (theta[1]+theta[2]) * es;
    resultdt.col(3) = theta[2] * es;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVaDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);


    resultDx.slice(0).col(0) = -theta[0] * s;
    resultDx.slice(0).col(1) = -theta[0] * e;
    resultDx.slice(0).col(2).fill(theta[1] + theta[2]);

    resultDx.slice(1).col(0) = -2*theta[0] * s;
    resultDx.slice(1).col(1) = -2*theta[0] * e;
    resultDx.slice(1).col(2).fill(2*theta[1]);

    resultDx.slice(2).col(0) = theta[0] * s;
    resultDx.slice(2).col(1) = theta[0] * e;
    resultDx.slice(2).col(2).fill(-theta[1] - theta[2]);

    resultDx.slice(3).col(2).fill(theta[2]);

    return resultDx;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVaDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & p = x.col(3);


    resultDtheta.slice(0).col(0) = -e % s;
    resultDtheta.slice(0).col(1) = es;
    resultDtheta.slice(0).col(2) = es;

    resultDtheta.slice(1).col(0) = -2 * e % s;
    resultDtheta.slice(1).col(1) = 2 * es;

    resultDtheta.slice(2).col(0) = e % s;
    resultDtheta.slice(2).col(1) = -es;
    resultDtheta.slice(2).col(2) = -es;

    resultDtheta.slice(3).col(2) = es;

    return resultDtheta;
}

// [[Rcpp::export]]
arma::mat MichaelisMentenModelVb6pODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(2);
    const double km2 = theta(3);
    const double k3 = theta(4);
    const double km3 = theta(5);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -k1 * e % s + km1 * es + k3 * es2 - km3 * e % p;
    resultdt.col(1) = -k1 * e % s - k2 * es % s + km1 * es;
    resultdt.col(2) = k1 * e % s - km1 * es - k2 * es % s + km2 * es2;
    resultdt.col(3) = k2 * es % s - (km2 + k3) * es2 + km3 * e %p;
    resultdt.col(4) = k3 * es2 - km3 * e % p;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVb6pDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(2);
    const double km2 = theta(3);
    const double k3 = theta(4);
    const double km3 = theta(5);

    resultDx.slice(0).col(0) = -k1 * s - km3 * p;
    resultDx.slice(0).col(1) = -k1 * e;
    resultDx.slice(0).col(2).fill(km1);
    resultDx.slice(0).col(3).fill(k3);
    resultDx.slice(0).col(4) = -km3 * e;

    resultDx.slice(1).col(0) = -k1 * s;
    resultDx.slice(1).col(1) = -k1*e - k2*es;
    resultDx.slice(1).col(2) = -k2*s + km1;

    resultDx.slice(2).col(0) = k1*s;
    resultDx.slice(2).col(1) = k1*e-k2*es;
    resultDx.slice(2).col(2) = -km1 - k2*s;
    resultDx.slice(2).col(3).fill(km2);

    resultDx.slice(3).col(0) = km3*p;
    resultDx.slice(3).col(1) = k2*es;
    resultDx.slice(3).col(2) = k2*s;
    resultDx.slice(3).col(3).fill(-(km2 + k3));
    resultDx.slice(3).col(4) = km3 * e;

    resultDx.slice(4).col(0) = -km3*p;
    resultDx.slice(4).col(3).fill(k3);
    resultDx.slice(4).col(4) = -km3 * e;

    return resultDx;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVb6pDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(2);
    const double km2 = theta(3);
    const double k3 = theta(4);
    const double km3 = theta(5);

    resultDtheta.slice(0).col(0) = -e % s;
    resultDtheta.slice(0).col(1) = es;
    resultDtheta.slice(0).col(4) = es2;
    resultDtheta.slice(0).col(5) = -e % p;

    resultDtheta.slice(1).col(0) = -e % s;
    resultDtheta.slice(1).col(1) = es;
    resultDtheta.slice(1).col(2) = -es % s;

    resultDtheta.slice(2).col(0) = e % s;
    resultDtheta.slice(2).col(1) = -es;
    resultDtheta.slice(2).col(2) = -es % s;
    resultDtheta.slice(2).col(3) = es2;

    resultDtheta.slice(3).col(2) = es % s;
    resultDtheta.slice(3).col(3) = -es2;
    resultDtheta.slice(3).col(4) = -es2;
    resultDtheta.slice(3).col(5) = e % p;

    resultDtheta.slice(4).col(4) = es2;
    resultDtheta.slice(4).col(5) = -e % p;

    return resultDtheta;
}
