// Fit ODE using Gaussian process
// simple case of FitzHugh-Nagumo (FN) Model
// assume all variables are observed with noise
// Fixed covar function: eta_sq=1, rho_sq=1, sigma=0.1

data {
  int<lower=1> N;
  real robs[N];
  real vobs[N];
  real drobs[N];
  real dvobs[N];
  real time[N];
}
transformed data {
  vector[N] mu;
  real delta = 1e-9;
  for (i in 1:N)
    mu[i] = 0;
  print("Finished transforming data.");
}
parameters {
  // vector[N1+N2+N3] f;
  // real<lower=0> abc[3];
  real<lower=0> sigma;
  real<lower=0> gamma;
  real<lower=0> rphi[2];
  real<lower=0> vphi[2];
  vector[N] reta;
  vector[N] veta;
}
model {
  matrix[N,N] C_rphi;
  matrix[N,N] L_C_rphi;
  matrix[N,N] inv_L_C_rphi;
  vector[N] rtrue;
  matrix[N,N] C_vphi;
  matrix[N,N] L_C_vphi;
  matrix[N,N] inv_L_C_vphi;
  vector[N] vtrue;
  
  vector[N] drtrue;
  vector[N] dvtrue;
  
  matrix[N,N] K_rphi;
  matrix[N,N] K_vphi;
  matrix[N,N] L_K_rphi;
  matrix[N,N] L_K_vphi;
  
  matrix[N,N] dC_rphi;
  matrix[N,N] ddC_rphi;
  matrix[N,N] dC_vphi;
  matrix[N,N] ddC_vphi;
  
  vector[N] m_rphi_rtrue;
  vector[N] m_vphi_vtrue;
  
  vector[N] dreta;
  vector[N] dveta;
  
  real r;
  real r2;
  
  for (i in 1:N)
    for (j in 1:N){
      r = fabs(time[i] - time[j]);
      r2 = pow(r, 2);
      C_rphi[i,j] = rphi[1] * (1 + ((sqrt(5)*r)/rphi[2]) + ((5*r2)/(3*pow(rphi[2],2)))) * exp((-sqrt(5)*r)/rphi[2]);
      C_vphi[i,j] = vphi[1] * (1 + ((sqrt(5)*r)/vphi[2]) + ((5*r2)/(3*pow(vphi[2],2)))) * exp((-sqrt(5)*r)/vphi[2]);
      if(i==j){
        C_rphi[i,j] = C_rphi[i,j] + delta;
        C_vphi[i,j] = C_vphi[i,j] + delta;
      }
      dC_rphi[i,j] = (2*step(j-i) - 1) * (rphi[1] * exp((-sqrt(5)*r)/rphi[2])) * 
        (((5*r)/(3*pow(rphi[2],2))) + ((5*sqrt(5)*r2)/(3*pow(rphi[2],3))));
      
      dC_vphi[i,j] = (2*step(j-i) - 1) * (vphi[1] * exp((-sqrt(5)*r)/vphi[2])) * 
        (((5*r)/(3*pow(vphi[2],2))) + ((5*sqrt(5)*r2)/(3*pow(vphi[2],3))));
      
      ddC_rphi[i,j] = (rphi[1]*exp((-sqrt(5)*r)/rphi[2])) * ((5/(3*pow(rphi[2],2))) 
        + ((5*sqrt(5)*r)/(3*pow(rphi[2],3))) - ((25*r2)/(3*pow(rphi[2],4))));
      ddC_vphi[i,j] = (vphi[1]*exp((-sqrt(5)*r)/vphi[2])) * ((5/(3*pow(vphi[2],2))) 
        + ((5*sqrt(5)*r)/(3*pow(vphi[2],3))) - ((25*r2)/(3*pow(vphi[2],4))));
    }
  
  L_C_rphi = cholesky_decompose(C_rphi);
  rtrue = L_C_rphi * reta;
  inv_L_C_rphi = inverse(L_C_rphi);
  
  L_C_vphi = cholesky_decompose(C_vphi);
  vtrue = L_C_vphi * veta;
  inv_L_C_vphi = inverse(L_C_vphi);
  
  m_rphi_rtrue = dC_rphi' * inv_L_C_rphi' * inv_L_C_rphi * rtrue;
  m_vphi_vtrue = dC_vphi' * inv_L_C_vphi' * inv_L_C_vphi * vtrue;
  
  K_rphi = inv_L_C_rphi * dC_rphi;
  K_rphi = ddC_rphi - K_rphi' * K_rphi;
  K_vphi = inv_L_C_vphi * dC_vphi;
  K_vphi = ddC_vphi - K_vphi' * K_vphi;
  for(i in 1:N){
    K_rphi[i,i] = K_rphi[i,i]+delta;
    K_vphi[i,i] = K_vphi[i,i]+delta;
  }
  L_K_rphi = cholesky_decompose(K_rphi);
  L_K_vphi = cholesky_decompose(K_vphi);
  
  drtrue = m_rphi_rtrue + K_rphi * dreta;
  dvtrue = m_vphi_vtrue + K_vphi * dveta;
  
  rphi[1] ~ cauchy(0,5);
  rphi[2] ~ cauchy(0,5);
  vphi[1] ~ cauchy(0,5);
  vphi[2] ~ cauchy(0,5);
  
  sigma ~ cauchy(0,5);
  gamma ~ cauchy(0,5);
  // abc ~ cauchy(0,5);
  
  reta ~ normal(0, 1);
  veta ~ normal(0, 1);
  
  dreta ~ normal(0, 1);
  dveta ~ normal(0, 1);
  
  robs ~ normal(rtrue, sigma);
  vobs ~ normal(vtrue, sigma);
  
  drobs ~ normal(drtrue, gamma);
  dvobs ~ normal(dvtrue, gamma);
}
