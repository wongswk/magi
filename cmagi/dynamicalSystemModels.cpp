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
arma::mat MichaelisMentenlog1xModelODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = arma::exp(x.col(0)) - 1;
    const vec & s = arma::exp(x.col(1)) - 1;
    const vec & es = arma::exp(x.col(2)) - 1;
    const vec & p = arma::exp(x.col(3)) - 1;

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = (-theta[0] * e % s + (theta[1]+theta[2]) * es) / (e + 1);
    resultdt.col(1) = (-theta[0] * e % s + (theta[1]) * es) / (s + 1);
    resultdt.col(2) = (theta[0] * e % s - (theta[1]+theta[2]) * es) / (es + 1);
    resultdt.col(3) = (theta[2] * es) / (p + 1);

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenlog1xModelDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = arma::exp(x.col(0)) - 1;
    const vec & s = arma::exp(x.col(1)) - 1;
    const vec & es = arma::exp(x.col(2)) - 1;
    const vec & p = arma::exp(x.col(3)) - 1;

    resultDx.slice(0).col(0) = (theta[0] * (e - 1) % s - (theta[1] + theta[2]) * es) / (e + 1);
    resultDx.slice(0).col(1) = -theta[0] * e % (s + 1) / (e + 1);
    resultDx.slice(0).col(2) = (theta[1] + theta[2]) * (es + 1) / (e + 1);

    resultDx.slice(1).col(0) = -theta[0] * s / (s + 1) % (e + 1);
    resultDx.slice(1).col(1) = -theta[0] * e - (-theta[0] * e % s + (theta[1]) * es) / (s + 1);
    resultDx.slice(1).col(2) = theta[1] / (s + 1) % (es + 1);

    resultDx.slice(2).col(0) = theta[0] * s / (es + 1) % (e + 1);
    resultDx.slice(2).col(1) = theta[0] * e / (es + 1) % (s + 1);
    resultDx.slice(2).col(2) = (-theta[1] - theta[2]) - (theta[0] * e % s - (theta[1]+theta[2]) * es) / (es + 1);

    resultDx.slice(3).col(2) = theta[2] / (p + 1) % (es + 1);
    resultDx.slice(3).col(3) = -(theta[2] * es) / (p + 1);

    return resultDx;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenlog1xModelDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & e = arma::exp(x.col(0)) - 1;
    const vec & s = arma::exp(x.col(1)) - 1;
    const vec & es = arma::exp(x.col(2)) - 1;
    const vec & p = arma::exp(x.col(3)) - 1;

    resultDtheta.slice(0).col(0) = (-e % s) / (e + 1);
    resultDtheta.slice(0).col(1) = es / (e + 1);
    resultDtheta.slice(0).col(2) = es / (e + 1);

    resultDtheta.slice(1).col(0) = (-e % s) / (s + 1);
    resultDtheta.slice(1).col(1) = es / (s + 1);

    resultDtheta.slice(2).col(0) = e % s / (es + 1);
    resultDtheta.slice(2).col(1) = -es / (es + 1);
    resultDtheta.slice(2).col(2) = -es / (es + 1);

    resultDtheta.slice(3).col(2) = es / (p + 1);

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


// [[Rcpp::export]]
arma::mat MichaelisMentenModelVb4pODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -k1 * e % s + km1 * es + k3 * es2 - km3 * e % p;
    resultdt.col(1) = -k1 * e % s - k2 * es % s + km1 * es;
    resultdt.col(2) = k1 * e % s - km1 * es - k2 * es % s + km2 * es2;
    resultdt.col(3) = k2 * es % s - (km2 + k3) * es2 + km3 * e %p;
    resultdt.col(4) = k3 * es2 - km3 * e % p;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVb4pDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);

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
arma::cube MichaelisMentenModelVb4pDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, 6, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);

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

    cube resultDtheta4p(x.n_rows, 4, x.n_cols, fill::zeros);

    for (int it = 0; it < 5; it++){
        resultDtheta4p.slice(it).col(0) = resultDtheta.slice(it).col(0) + resultDtheta.slice(it).col(2);
        resultDtheta4p.slice(it).col(1) = resultDtheta.slice(it).col(1) + resultDtheta.slice(it).col(3);
        resultDtheta4p.slice(it).col(2) = resultDtheta.slice(it).col(4);
        resultDtheta4p.slice(it).col(3) = resultDtheta.slice(it).col(5);
    }
    
    return resultDtheta4p;
}

// [[Rcpp::export]]
arma::mat MichaelisMentenModelVb2pODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = 0;
    const double k2 = theta(0);
    const double km2 = 0;
    const double k3 = theta(1);
    const double km3 = 0;

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -k1 * e % s + km1 * es + k3 * es2 - km3 * e % p;
    resultdt.col(1) = -k1 * e % s - k2 * es % s + km1 * es;
    resultdt.col(2) = k1 * e % s - km1 * es - k2 * es % s + km2 * es2;
    resultdt.col(3) = k2 * es % s - (km2 + k3) * es2 + km3 * e %p;
    resultdt.col(4) = k3 * es2 - km3 * e % p;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVb2pDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = 0;
    const double k2 = theta(0);
    const double km2 = 0;
    const double k3 = theta(1);
    const double km3 = 0;

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
arma::cube MichaelisMentenModelVb2pDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, 6, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = 0;
    const double k2 = theta(0);
    const double km2 = 0;
    const double k3 = theta(1);
    const double km3 = 0;

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

    cube resultDtheta2p(x.n_rows, 2, x.n_cols, fill::zeros);

    for (int it = 0; it < 5; it++){
        resultDtheta2p.slice(it).col(0) = resultDtheta.slice(it).col(0) + resultDtheta.slice(it).col(2);
        resultDtheta2p.slice(it).col(1) = resultDtheta.slice(it).col(4);
    }

    return resultDtheta2p;
}

// [[Rcpp::export]]
arma::mat MichaelisMentenModelVc7pODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
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
    const double kin = theta(6);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -(k1+kin) * e % s + km1 * es + k3 * es2 - km3 * e % p;
    resultdt.col(1) = -k1 * e % s - k2 * es % s + km1 * es;
    resultdt.col(2) = k1 * e % s - km1 * es - k2 * es % s + km2 * es2;
    resultdt.col(3) = k2 * es % s - (km2 + k3) * es2 + km3 * e %p;
    resultdt.col(4) = k3 * es2 - km3 * e % p;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVc7pDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
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
    const double kin = theta(6);

    resultDx.slice(0).col(0) = -(k1 + kin) * s - km3 * p;
    resultDx.slice(0).col(1) = -(k1 + kin) * e;
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
arma::cube MichaelisMentenModelVc7pDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
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
    const double kin = theta(6);

    resultDtheta.slice(0).col(0) = -e % s;
    resultDtheta.slice(0).col(1) = es;
    resultDtheta.slice(0).col(4) = es2;
    resultDtheta.slice(0).col(5) = -e % p;
    resultDtheta.slice(0).col(6) = -e % s;

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


// [[Rcpp::export]]
arma::mat MichaelisMentenModelVc5pODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);
    const double kin = theta(4);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -(k1+kin) * e % s + km1 * es + k3 * es2 - km3 * e % p;
    resultdt.col(1) = -k1 * e % s - k2 * es % s + km1 * es;
    resultdt.col(2) = k1 * e % s - km1 * es - k2 * es % s + km2 * es2;
    resultdt.col(3) = k2 * es % s - (km2 + k3) * es2 + km3 * e %p;
    resultdt.col(4) = k3 * es2 - km3 * e % p;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube MichaelisMentenModelVc5pDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);
    const double kin = theta(4);

    resultDx.slice(0).col(0) = -(k1 + kin) * s - km3 * p;
    resultDx.slice(0).col(1) = -(k1 + kin) * e;
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
arma::cube MichaelisMentenModelVc5pDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, 7, x.n_cols, fill::zeros);

    const vec & e = x.col(0);
    const vec & s = x.col(1);
    const vec & es = x.col(2);
    const vec & es2 = x.col(3);
    const vec & p = x.col(4);

    const double k1 = theta(0);
    const double km1 = theta(1);
    const double k2 = theta(0);
    const double km2 = theta(1);
    const double k3 = theta(2);
    const double km3 = theta(3);
    const double kin = theta(4);

    resultDtheta.slice(0).col(0) = -e % s;
    resultDtheta.slice(0).col(1) = es;
    resultDtheta.slice(0).col(4) = es2;
    resultDtheta.slice(0).col(5) = -e % p;
    resultDtheta.slice(0).col(6) = -e % s;

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

    cube resultDtheta5p(x.n_rows, 5, x.n_cols, fill::zeros);

    for (int it = 0; it < 5; it++){
        resultDtheta5p.slice(it).col(0) = resultDtheta.slice(it).col(0) + resultDtheta.slice(it).col(2);
        resultDtheta5p.slice(it).col(1) = resultDtheta.slice(it).col(1) + resultDtheta.slice(it).col(3);
        resultDtheta5p.slice(it).col(2) = resultDtheta.slice(it).col(4);
        resultDtheta5p.slice(it).col(3) = resultDtheta.slice(it).col(5);
        resultDtheta5p.slice(it).col(4) = resultDtheta.slice(it).col(6);
    }

    return resultDtheta5p;
}

// [[Rcpp::export]]
arma::mat lacOperonODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & ri = x.col(0);
    const vec & i = x.col(1);
    const vec & lactose = x.col(2);
    const vec & ilactose = x.col(3);
    const vec & op = x.col(4);
    const vec & iop = x.col(5);
    const vec & rnap = x.col(6);
    const vec & rnapo = x.col(7);
    const vec & r = x.col(8);
    const vec & z = x.col(9);

//    hard code [i] component to 1
//    const double iconstant = theta(0);
    const double iconstant = 1.0;

    const vec & k = theta;

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = k(1) * iconstant - k(12) * ri;
    resultdt.col(1) = k(2) * ri - k(3) * i % lactose + k(4) * ilactose - k(5) * i % op + k(6) * iop - k(13) * i;
    resultdt.col(2) = k(4) * ilactose - k(3) * i % lactose + k(14) * ilactose - k(11) * lactose % z;
    resultdt.col(3) = k(3) * i % lactose - k(4) * ilactose - k(14) * ilactose;
    resultdt.col(4) = k(6) * iop - k(5) * i % op - k(7) * op % rnap + (k(8) + k(9)) * rnapo;
    resultdt.col(5) = k(5) * i % op - k(6) * iop;
    resultdt.col(6) = (k(8) + k(9)) * rnapo - k(7) * op % rnap;
    resultdt.col(7) = k(7) * op % rnap - (k(8)+k(9)) * rnapo;
    resultdt.col(8) = k(9) * rnapo - k(15) * r;
    resultdt.col(9) = k(10)*r - k(16) * z;

    return resultdt;
}

// [[Rcpp::export]]
arma::cube lacOperonDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & ri = x.col(0);
    const vec & i = x.col(1);
    const vec & lactose = x.col(2);
    const vec & ilactose = x.col(3);
    const vec & op = x.col(4);
    const vec & iop = x.col(5);
    const vec & rnap = x.col(6);
    const vec & rnapo = x.col(7);
    const vec & r = x.col(8);
    const vec & z = x.col(9);

//    hard code [i] component to 1
//    const double iconstant = theta(0);
    const double iconstant = 1.0;
    const vec & k = theta;

    resultDx.slice(0).col(0).fill(-k(12));

    resultDx.slice(1).col(0).fill(k(2));
    resultDx.slice(1).col(1) = -k(3) * lactose - k(5) * op - k(13);
    resultDx.slice(1).col(2) = -k(3) * i;
    resultDx.slice(1).col(3).fill(k(4));
    resultDx.slice(1).col(4) = -k(5) * i;
    resultDx.slice(1).col(5).fill(k(6));

    resultDx.slice(2).col(1) = -k(3)*lactose;
    resultDx.slice(2).col(2) = -k(3)*i - k(11)*z;
    resultDx.slice(2).col(3).fill(k(4) + k(14));
    resultDx.slice(2).col(9) = -k(11)*lactose;

    resultDx.slice(3).col(1) = k(3)*lactose;
    resultDx.slice(3).col(2) = k(3)*i;
    resultDx.slice(3).col(3).fill(-k(4) - k(14));

    resultDx.slice(4).col(1) = -k(5) * op;
    resultDx.slice(4).col(4) = -k(5) * i - k(7)*rnap;
    resultDx.slice(4).col(5).fill(k(6));
    resultDx.slice(4).col(6) = -k(7)*op;
    resultDx.slice(4).col(7).fill(k(8) + k(9));

    resultDx.slice(5).col(1) = k(5) * op;
    resultDx.slice(5).col(4) = k(5) * i;
    resultDx.slice(5).col(5).fill(-k(6));

    resultDx.slice(6).col(4) = -k(7)*rnap;
    resultDx.slice(6).col(6) = -k(7)*op;
    resultDx.slice(6).col(7).fill(k(8) + k(9));

    resultDx.slice(7).col(4) = k(7) * rnap;
    resultDx.slice(7).col(6) = k(7) * op;
    resultDx.slice(7).col(7).fill(-(k(8) + k(9)));

    resultDx.slice(8).col(7).fill(k(9));
    resultDx.slice(8).col(8).fill(-k(15));

    resultDx.slice(9).col(8).fill(k(10));
    resultDx.slice(9).col(9).fill(-k(16));

    return resultDx;
}

// [[Rcpp::export]]
arma::cube lacOperonDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & ri = x.col(0);
    const vec & i = x.col(1);
    const vec & lactose = x.col(2);
    const vec & ilactose = x.col(3);
    const vec & op = x.col(4);
    const vec & iop = x.col(5);
    const vec & rnap = x.col(6);
    const vec & rnapo = x.col(7);
    const vec & r = x.col(8);
    const vec & z = x.col(9);

//    hard code [i] component to 1
//    const double iconstant = theta(0);
    const double iconstant = 1.0;
    const vec & k = theta;

//    hard code [i] component to 1
//    resultDtheta.slice(0).col(0).fill(k(1));
    resultDtheta.slice(0).col(0).fill(0);
    resultDtheta.slice(0).col(1).fill(iconstant);
    resultDtheta.slice(0).col(12) = -ri;

    resultDtheta.slice(1).col(2) = ri;
    resultDtheta.slice(1).col(3) = -i%lactose;
    resultDtheta.slice(1).col(4) = ilactose;
    resultDtheta.slice(1).col(5) = - i%op;
    resultDtheta.slice(1).col(6) = iop;
    resultDtheta.slice(1).col(13) = -i;

    resultDtheta.slice(2).col(4) = ilactose;
    resultDtheta.slice(2).col(3) = -i%lactose;
    resultDtheta.slice(2).col(14) = ilactose;
    resultDtheta.slice(2).col(11) = -lactose%z;

    resultDtheta.slice(3).col(3) = i%lactose;
    resultDtheta.slice(3).col(4) = -ilactose;
    resultDtheta.slice(3).col(14) = -ilactose;

    resultDtheta.slice(4).col(6) = iop;
    resultDtheta.slice(4).col(5) = -i%op;
    resultDtheta.slice(4).col(7) = -op%rnap;
    resultDtheta.slice(4).col(8) = rnapo;
    resultDtheta.slice(4).col(9) = rnapo;

    resultDtheta.slice(5).col(5) = i%op;
    resultDtheta.slice(5).col(6) = -iop;

    resultDtheta.slice(6).col(8) = rnapo;
    resultDtheta.slice(6).col(9) = rnapo;
    resultDtheta.slice(6).col(7) = -op%rnap;

    resultDtheta.slice(7).col(7) = op%rnap;
    resultDtheta.slice(7).col(8) = -rnapo;
    resultDtheta.slice(7).col(9) = -rnapo;

    resultDtheta.slice(8).col(9) = rnapo;
    resultDtheta.slice(8).col(15) = -r;

    resultDtheta.slice(9).col(10) = r;
    resultDtheta.slice(9).col(16) = -z;

    return resultDtheta;
}

// [[Rcpp::export]]
arma::mat repressilatorGeneRegulationODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & m_laci = x.col(0);
    const vec & m_tetr = x.col(1);
    const vec & m_ci = x.col(2);
    const vec & p_laci = x.col(3);
    const vec & p_tetr = x.col(4);
    const vec & p_ci = x.col(5);

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -m_laci + alpha / (1 + arma::pow(p_ci, n)) + alpha0;
    resultdt.col(1) = -m_tetr + alpha / (1 + arma::pow(p_laci, n)) + alpha0;
    resultdt.col(2) = -m_ci + alpha / (1 + arma::pow(p_tetr, n)) + alpha0;
    resultdt.col(3) = -beta*(p_laci - m_laci);
    resultdt.col(4) = -beta*(p_tetr - m_tetr);
    resultdt.col(5) = -beta*(p_ci - m_ci);

    return resultdt;
}


// [[Rcpp::export]]
arma::cube repressilatorGeneRegulationDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & m_laci = x.col(0);
    const vec & m_tetr = x.col(1);
    const vec & m_ci = x.col(2);
    const vec & p_laci = x.col(3);
    const vec & p_tetr = x.col(4);
    const vec & p_ci = x.col(5);

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    resultDx.slice(0).col(0).fill(-1);
    resultDx.slice(0).col(5) = alpha * (-1) * arma::pow(1 + arma::pow(p_ci, n), -2) * n % arma::pow(p_ci, n-1);
    resultDx.slice(1).col(1).fill(-1);
    resultDx.slice(1).col(3) = alpha * (-1) * arma::pow(1 + arma::pow(p_laci, n), -2) * n % arma::pow(p_laci, n-1);
    resultDx.slice(2).col(2).fill(-1);
    resultDx.slice(2).col(4) = alpha * (-1) * arma::pow(1 + arma::pow(p_tetr, n), -2) * n % arma::pow(p_tetr, n-1);

    resultDx.slice(3).col(0).fill(beta);
    resultDx.slice(3).col(3).fill(-beta);
    resultDx.slice(4).col(1).fill(beta);
    resultDx.slice(4).col(4).fill(-beta);
    resultDx.slice(5).col(2).fill(beta);
    resultDx.slice(5).col(5).fill(-beta);

    return resultDx;
}


// [[Rcpp::export]]
arma::cube repressilatorGeneRegulationDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & m_laci = x.col(0);
    const vec & m_tetr = x.col(1);
    const vec & m_ci = x.col(2);
    const vec & p_laci = x.col(3);
    const vec & p_tetr = x.col(4);
    const vec & p_ci = x.col(5);

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    resultDtheta.slice(0).col(0).fill(1);
    resultDtheta.slice(0).col(1) = 1 / (1 + arma::pow(p_ci, n));
    resultDtheta.slice(0).col(2) = alpha * (-1) * pow(1 + arma::pow(p_ci, n), -2) % arma::pow(p_ci, n) % arma::log(p_ci);
    resultDtheta.slice(1).col(0).fill(1);
    resultDtheta.slice(1).col(1) = 1 / (1 + arma::pow(p_laci, n));
    resultDtheta.slice(1).col(2) = alpha * (-1) * pow(1 + arma::pow(p_laci, n), -2) % arma::pow(p_laci, n) % arma::log(p_laci);
    resultDtheta.slice(2).col(0).fill(1);
    resultDtheta.slice(2).col(1) = 1 / (1 + arma::pow(p_tetr, n));
    resultDtheta.slice(2).col(2) = alpha * (-1) * pow(1 + arma::pow(p_tetr, n), -2) % arma::pow(p_tetr, n) % arma::log(p_tetr);

    resultDtheta.slice(3).col(3) = x.col(0) - x.col(3);
    resultDtheta.slice(4).col(3) = x.col(1) - x.col(4);
    resultDtheta.slice(5).col(3) = x.col(2) - x.col(5);

    return resultDtheta;
}

// [[Rcpp::export]]
arma::mat repressilatorGeneRegulationLogODE(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    const vec & m_laci = arma::exp(x.col(0));
    const vec & m_tetr = arma::exp(x.col(1));
    const vec & m_ci = arma::exp(x.col(2));
    const vec & p_laci = arma::exp(x.col(3));
    const vec & p_tetr = arma::exp(x.col(4));
    const vec & p_ci = arma::exp(x.col(5));

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    mat resultdt(x.n_rows, x.n_cols);

    resultdt.col(0) = -1 + (alpha / (1 + arma::pow(p_ci, n)) + alpha0) / m_laci;
    resultdt.col(1) = -1 + (alpha / (1 + arma::pow(p_laci, n)) + alpha0) / m_tetr;
    resultdt.col(2) = -1 + (alpha / (1 + arma::pow(p_tetr, n)) + alpha0) / m_ci;
    resultdt.col(3) = (-beta*(1 - m_laci / p_laci));
    resultdt.col(4) = (-beta*(1 - m_tetr / p_tetr)) ;
    resultdt.col(5) = (-beta*(1 - m_ci / p_ci));

    return resultdt;
}


// [[Rcpp::export]]
arma::cube repressilatorGeneRegulationLogDx(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDx(x.n_rows, x.n_cols, x.n_cols, fill::zeros);

    const vec & tm_laci = x.col(0);
    const vec & tm_tetr = x.col(1);
    const vec & tm_ci = x.col(2);
    const vec & tp_laci = x.col(3);
    const vec & tp_tetr = x.col(4);
    const vec & tp_ci = x.col(5);

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    resultDx.slice(0).col(0) = -(alpha / (1 + arma::exp(n * tp_ci)) + alpha0) % arma::exp(-tm_laci);
    resultDx.slice(0).col(5) = alpha * arma::exp(-tm_laci) * (-1) % arma::pow(1 + arma::exp(n * tp_ci), -2) * n % arma::exp(n * tp_ci);
    resultDx.slice(1).col(1) = -(alpha / (1 + arma::exp(n * tp_laci)) + alpha0) % arma::exp(-tm_tetr);
    resultDx.slice(1).col(3) = alpha * arma::exp(-tm_tetr) * (-1) % arma::pow(1 + arma::exp(n * tp_laci), -2) * n % arma::exp(n * tp_laci);
    resultDx.slice(2).col(2) = -(alpha / (1 + arma::exp(n * tp_tetr)) + alpha0) % arma::exp(-tm_ci);
    resultDx.slice(2).col(4) = alpha * arma::exp(-tm_ci) * (-1) % arma::pow(1 + arma::exp(n * tp_tetr), -2) * n % arma::exp(n * tp_tetr);

    resultDx.slice(3).col(0) = beta * arma::exp(tm_laci - tp_laci);
    resultDx.slice(3).col(3) = -beta * arma::exp(tm_laci - tp_laci);;
    resultDx.slice(4).col(1) = beta * arma::exp(tm_tetr - tp_tetr);
    resultDx.slice(4).col(4) = -beta * arma::exp(tm_tetr - tp_tetr);
    resultDx.slice(5).col(2) = beta * arma::exp(tm_ci - tp_ci);
    resultDx.slice(5).col(5) = -beta * arma::exp(tm_ci - tp_ci);

    return resultDx;
}


// [[Rcpp::export]]
arma::cube repressilatorGeneRegulationLogDtheta(const arma::vec & theta, const arma::mat & x, const arma::vec & tvec) {
    cube resultDtheta(x.n_rows, theta.size(), x.n_cols, fill::zeros);

    const vec & tm_laci = x.col(0);
    const vec & tm_tetr = x.col(1);
    const vec & tm_ci = x.col(2);
    const vec & tp_laci = x.col(3);
    const vec & tp_tetr = x.col(4);
    const vec & tp_ci = x.col(5);

    const vec p_ci = arma::exp(tp_ci);
    const vec p_laci = arma::exp(tp_laci);
    const vec p_tetr = arma::exp(tp_tetr);

    const double alpha0 = theta(0);
    const double alpha = theta(1);
    const double n = theta(2);
    const double beta = theta(3);

    resultDtheta.slice(0).col(0) = arma::exp(-x.col(0));
    resultDtheta.slice(0).col(1) = 1 / (1 + arma::exp(n * tp_ci)) % arma::exp(-x.col(0));
    resultDtheta.slice(0).col(2) = alpha * arma::exp(-x.col(0)) * (-1) % pow(1 + arma::pow(p_ci, n), -2) % arma::pow(p_ci, n) % arma::log(p_ci);
    resultDtheta.slice(1).col(0) = arma::exp(-x.col(1));
    resultDtheta.slice(1).col(1) = 1 / (1 + arma::exp(n * tp_laci)) % arma::exp(-x.col(1));
    resultDtheta.slice(1).col(2) = alpha * arma::exp(-x.col(1)) * (-1) % pow(1 + arma::pow(p_laci, n), -2) % arma::pow(p_laci, n) % arma::log(p_laci);
    resultDtheta.slice(2).col(0) = arma::exp(-x.col(2));
    resultDtheta.slice(2).col(1) = 1 / (1 + arma::exp(n * tp_tetr)) % arma::exp(-x.col(2));
    resultDtheta.slice(2).col(2) = alpha * arma::exp(-x.col(2)) * (-1) % pow(1 + arma::pow(p_tetr, n), -2) % arma::pow(p_tetr, n) % arma::log(p_tetr);

    resultDtheta.slice(3).col(3) = arma::exp(x.col(0) - x.col(3)) - 1;
    resultDtheta.slice(4).col(3) = arma::exp(x.col(1) - x.col(4)) - 1;
    resultDtheta.slice(5).col(3) = arma::exp(x.col(2) - x.col(5)) - 1;

    return resultDtheta;
}
