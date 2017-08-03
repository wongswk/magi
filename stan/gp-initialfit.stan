// Fit ODE using Gaussian process
// simple case of FitzHugh-Nagumo (FN) Model
// assume all variables are observed with noise
// Fixed covar function: eta_sq=1, rho_sq=1, sigma=0.1

data {
  int<lower=1> N;
  real robs[N];
  real vobs[N];
  real time[N];
}
transformed data {
  vector[N] mu;
  // real gamma = 0;
  real delta = 1e-9;
  // real sigma = 0.1;
  matrix[N,N] r;
  matrix[N,N] r2;
  for (i in 1:N){
    mu[i] = 0;
    for (j in 1:N){
      r[i,j] = fabs(time[i] - time[j]);
      r2[i,j] = pow(r[i,j], 2);
    }
  }
}
parameters {
  // vector[N1+N2+N3] f;
  // real<lower=0> abc[3];
  real<lower=0> sigma;
  real<lower=0> rphi[2];
  real<lower=0> vphi[2];
  vector[N] reta;
  vector[N] veta;
}
model {
  matrix[N,N] C_rphi;
  matrix[N,N] L_C_rphi;
  vector[N] rtrue;
  matrix[N,N] C_vphi;
  matrix[N,N] L_C_vphi;
  vector[N] vtrue;
  for (i in 1:N)
    for (j in 1:N){
      C_rphi[i,j] = rphi[1] * (1 + ((sqrt(5)*r[i,j])/rphi[2]) + ((5*r2[i,j])/(3*pow(rphi[2],2)))) * 
        exp((-sqrt(5)*r[i,j])/rphi[2]);
      C_vphi[i,j] = vphi[1] * (1 + ((sqrt(5)*r[i,j])/vphi[2]) + ((5*r2[i,j])/(3*pow(vphi[2],2)))) * 
        exp((-sqrt(5)*r[i,j])/vphi[2]);
      if(i==j){
        C_rphi[i,j] = C_rphi[i,j] + delta;
        C_vphi[i,j] = C_vphi[i,j] + delta;
      }
    }
  
  L_C_rphi = cholesky_decompose(C_rphi);
  rtrue = L_C_rphi * reta;
  
  L_C_vphi = cholesky_decompose(C_vphi);
  vtrue = L_C_vphi * veta;
  
  rphi[1] ~ cauchy(0,5);
  rphi[2] ~ cauchy(0,5);
  vphi[1] ~ cauchy(0,5);
  vphi[2] ~ cauchy(0,5);
  
  sigma ~ cauchy(0,5);
  
  reta ~ normal(0, 1);
  veta ~ normal(0, 1);
  
  robs ~ normal(rtrue, sigma);
  vobs ~ normal(vtrue, sigma);
}

generated quantities {
  matrix[N,N] C_rphi;
  matrix[N,N] L_C_rphi;
  vector[N] rtrue;
  matrix[N,N] C_vphi;
  matrix[N,N] L_C_vphi;
  vector[N] vtrue;
  for (i in 1:N)
    for (j in 1:N){
      C_rphi[i,j] = rphi[1] * (1 + ((sqrt(5)*r[i,j])/rphi[2]) + ((5*r2[i,j])/(3*pow(rphi[2],2)))) * 
        exp((-sqrt(5)*r[i,j])/rphi[2]);
      C_vphi[i,j] = vphi[1] * (1 + ((sqrt(5)*r[i,j])/vphi[2]) + ((5*r2[i,j])/(3*pow(vphi[2],2)))) * 
        exp((-sqrt(5)*r[i,j])/vphi[2]);
      if(i==j){
        C_rphi[i,j] = C_rphi[i,j] + delta;
        C_vphi[i,j] = C_vphi[i,j] + delta;
      }
    }
  
  L_C_rphi = cholesky_decompose(C_rphi);
  rtrue = L_C_rphi * reta;
  
  L_C_vphi = cholesky_decompose(C_vphi);
  vtrue = L_C_vphi * veta;
}