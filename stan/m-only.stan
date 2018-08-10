// Fit ODE using Finite difference
// simple case of FitzHugh-Nagumo (FN) Model
// assume all variables are observed with noise
// assume sigma is known
// Fixed covar function: eta_sq=1, rho_sq=1, sigma=0.1

data {
  int<lower=1> n_discret;
  int<lower=1> n_obs;
  real robs[n_obs];
  real vobs[n_obs];
  real time[n_discret];
  int obs_index[n_obs];
  real sigma_obs;
  real sigma_xdot;
  matrix[n_discret, n_discret] mphi_r;
  matrix[n_discret, n_discret] mphi_v;
}

parameters {
  // vector[N1+N2+N3] f;
  real<lower=0> abc[3];
  vector[n_discret] rtrue;
  vector[n_discret] vtrue;
}

model {
  vector[n_discret] drobs;
  vector[n_discret] dvobs;
  
  vector[n_discret] m_rtrue;
  vector[n_discret] m_vtrue;
  
  vector[n_obs] rtrue_at_obs;
  vector[n_obs] vtrue_at_obs;
  
  for (i in 1:n_discret){
    dvobs[i] = abc[3] * (vtrue[i] - pow(vtrue[i],3)/3.0 + rtrue[i]);  
    drobs[i] = -1.0/abc[3] * (vtrue[i] - abc[1] + abc[2]*rtrue[i]);
  }
  
  m_rtrue = mphi_r * rtrue;
  m_vtrue = mphi_v * vtrue;
  
  for (i in 1:n_obs){
    rtrue_at_obs[i] = rtrue[obs_index[i]];
    vtrue_at_obs[i] = vtrue[obs_index[i]];
  }
  
  abc[1] ~ cauchy(0,5);
  abc[2] ~ cauchy(0,5);
  abc[3] ~ cauchy(0,5);
  
  robs ~ normal(rtrue_at_obs, sigma_obs);
  vobs ~ normal(vtrue_at_obs, sigma_obs);
  
  drobs ~ normal(m_rtrue, sigma_xdot);
  dvobs ~ normal(m_vtrue, sigma_xdot);
}
