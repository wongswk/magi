// Fit ODE using Gaussian process
// simple case of FitzHugh-Nagumo (FN) Model
// assume all variables are observed with noise
// Fixed covar function: eta_sq=1, rho_sq=1, sigma_sq=0.1

data {
  int<lower=1> N;
  real robs[N];
  real time[N];
}
transformed data {
  vector[N] mu;
  real delta = 1e-9;
  for (i in 1:N)
    mu[i] = 0;
}
parameters {
  // vector[N1+N2+N3] f;
  // real<lower=0> abc[3];
  real<lower=0> sigma_sq;
  real<lower=0> rphi[2];
  vector[N] eta;
}
model {
  matrix[N,N] C_rphi;
  matrix[N,N] L_C_rphi;
  vector[N] rtrue;
  real r;
  real r2;
  for (i in 1:N)
    for (j in 1:N){
      r = fabs(time[i] - time[j]);
      r2 = pow(r, 2);
      C_rphi[i,j] = rphi[1] * (1 + ((sqrt(5)*r)/rphi[2]) + ((5*r2)/(3*pow(rphi[2],2)))) * exp((-sqrt(5)*r)/rphi[2]);
      if(i==j){
        C_rphi[i,j] = C_rphi[i,j] + delta;
      }
    }
  
  L_C_rphi = cholesky_decompose(C_rphi);
  rtrue = L_C_rphi * eta;
  
  rphi[1] ~ cauchy(0,5);
  rphi[2] ~ cauchy(0,5);
  sigma_sq ~ cauchy(0,5);
  
  eta ~ normal(0, 1);
  robs ~ normal(rtrue, sigma_sq);
}
