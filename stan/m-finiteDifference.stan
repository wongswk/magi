// Fit ODE using Finite difference
// simple case of FitzHugh-Nagumo (FN) Model
// assume all variables are observed with noise
// assume sigma is known
// Fixed covar function: eta_sq=1, rho_sq=1, sigma=0.1

data {
  int<lower=1> N;
  real robs[N];
  real vobs[N];
  real time[N];
  real sigma_obs;
  real sigma_xdot;
}

parameters {
  // vector[N1+N2+N3] f;
  real<lower=0> abc[3];
  vector[N] rtrue;
  vector[N] vtrue;
}

model {
  vector[N] drobs;
  vector[N] dvobs;
  
  vector[N] m_rtrue;
  vector[N] m_vtrue;
  
  for (i in 1:N){
    dvobs[i] = abc[3] * (vtrue[i] - pow(vtrue[i],3)/3.0 + rtrue[i]);  
    drobs[i] = -1.0/abc[3] * (vtrue[i] - abc[1] + abc[2]*rtrue[i]);
    
    if(i == 1){
      m_vtrue[i] = (vtrue[2] - vtrue[1]) / (time[2] - time[1]);
      m_rtrue[i] = (rtrue[2] - rtrue[1]) / (time[2] - time[1]);
    }else if(i == N){
      m_vtrue[i] = (vtrue[N] - vtrue[N-1]) / (time[N] - time[N-1]);
      m_rtrue[i] = (rtrue[N] - rtrue[N-1]) / (time[N] - time[N-1]);
    }else{
      m_vtrue[i] = (vtrue[i+1] - vtrue[i-1]) / (time[i+1] - time[i-1]);
      m_rtrue[i] = (rtrue[i+1] - rtrue[i-1]) / (time[i+1] - time[i-1]);
    }
  }
  
  abc[1] ~ cauchy(0,5);
  abc[2] ~ cauchy(0,5);
  abc[3] ~ cauchy(0,5);
  
  robs ~ normal(rtrue, sigma_obs);
  vobs ~ normal(vtrue, sigma_obs);
  
  drobs ~ normal(m_rtrue, sigma_xdot);
  dvobs ~ normal(m_vtrue, sigma_xdot);
}
