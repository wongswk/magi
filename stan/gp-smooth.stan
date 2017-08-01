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
  for (i in 1:N)
    mu[i] = 0;
}
parameters {
  // vector[N1+N2+N3] f;
  vector[N] rtrue;
  // real<lower=0> abc[3];
  real<lower=0> sigma_sq;
  real<lower=0> rphi[2];
}
model {
  matrix[N,N] C_rphi;
  real r;
  real r2;
  for (i in 1:N)
    for (j in 1:N){
      r = fabs(time[i] - time[j]);
      r2 = pow(r, 2);
      C_rphi[i,j] = rphi[1] * (1 + ((sqrt(5)*r)/rphi[2]) + ((5*r2)/(3*pow(rphi[2],2)))) * exp((-sqrt(5)*r)/rphi[2]);
    }
    
  rphi[1] ~ cauchy(0,5);
  rphi[2] ~ cauchy(0,5);
  
  rtrue ~ multi_normal(mu, C_rphi);
  robs ~ normal(rtrue,sigma_sq);
}
