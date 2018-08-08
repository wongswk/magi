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
  int finiteDifferenceType;
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
    
    if(i == 1){
      m_vtrue[i] = (vtrue[2] - vtrue[1]) / (time[2] - time[1]);
      m_rtrue[i] = (rtrue[2] - rtrue[1]) / (time[2] - time[1]);
    }else if(i == n_discret){
      m_vtrue[i] = (vtrue[n_discret] - vtrue[n_discret-1]) / (time[n_discret] - time[n_discret-1]);
      m_rtrue[i] = (rtrue[n_discret] - rtrue[n_discret-1]) / (time[n_discret] - time[n_discret-1]);
    }else{
      if(finiteDifferenceType == 0){
        m_vtrue[i] = (vtrue[i+1] - vtrue[i-1]) / (time[i+1] - time[i-1]);
        m_rtrue[i] = (rtrue[i+1] - rtrue[i-1]) / (time[i+1] - time[i-1]);
      }else if(finiteDifferenceType == 1){
        m_vtrue[i] = (vtrue[i+1] - vtrue[i]) / (time[i+1] - time[i]);
        m_rtrue[i] = (rtrue[i+1] - rtrue[i]) / (time[i+1] - time[i]);
      }
    }
  }
  
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
