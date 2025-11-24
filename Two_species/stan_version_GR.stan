functions {
  vector ode_function(real t,
                      vector y,
                      real alpha_tau,
                      real alpha_omega,
                      // real nu_tau,
                      real nu_omega,
                      real beta_tau_omega,
                      real beta_omega_tau,
                      real kappa_tau,
                      real kappa_omega
  ) {
  vector[2] dydt; // dydt[1] = tagged strain. dydt[2] = untagged strain.
  dydt[1] = alpha_tau * (y[1]) * log(kappa_tau*kappa_omega/(kappa_omega*y[1]+beta_tau_omega*kappa_tau*y[2]));
  dydt[2] = alpha_omega * (y[2]) * (1 - (y[2] / kappa_omega + beta_omega_tau * (y[1] / kappa_tau)) ^ nu_omega);
    
  return(dydt);
  }
}



// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> T; // Number of time points (including zero)
  int<lower=1> n_wells; 
  array[2*T] real ts; // Vector of times
  array[T] int t_index_g; // index of where the output times for g match the model times vector
  array[T] int t_index_b;
  array[T,n_wells] vector[2] y; // 3D array of measurements by time, well, and type
  array[n_wells] vector[2] y0; // Initial conditions by well and type
  real rho_1;
  real rho_2;
  int t_star;
}

// The parameters accepted by the model.
parameters {
  real<lower=0> sigma;
  real<lower=270,upper=380> delta_mean;
  real<lower=0.088,upper = 0.11> epsilon_mean;
  real<lower=0> delta_sd;
  real<lower=0> epsilon_sd;
  real<lower=0> c_tau_mean;
  real<lower=0> c_omega_mean;
  
  real<lower=130, upper=240> rho_omega;

  real<lower=0> kappa_tau_mean;
  real<lower=0> kappa_omega_mean;
  real<lower=0> c_tau_sd;
  real<lower=0> c_omega_sd;

  real<lower=0> kappa_tau_sd;
  real<lower=0> kappa_omega_sd;
  real<lower=0> alpha_tau_mean;
  real<lower=0> alpha_omega_mean;
  real<lower=0> alpha_tau_sd;
  real<lower=0> alpha_omega_sd;
  
  real<lower=0, upper=1> nu_tau;
  real<lower=0, upper=1> nu_omega;
  real<lower=0> beta_tau_omega;
  real<lower=0> beta_omega_tau;
  vector<lower=0>[n_wells] c_tau;
  vector<lower=0>[n_wells] c_omega;

  vector<lower=0>[n_wells] kappa_tau;
  vector<lower=0>[n_wells] kappa_omega;
  
  vector<lower=0>[n_wells] alpha_tau;
  vector<lower=0>[n_wells] alpha_omega;
  
  vector<lower=270,upper=380>[n_wells] delta;
  vector<lower=0.088,upper = 0.11>[n_wells] epsilon;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // Define priors
  sigma ~ exponential(1000);
  delta_mean ~ normal(331.16,12.94);
  epsilon_mean ~ normal(0.09555,0.001723);
  delta_sd ~ exponential(1);
  epsilon_sd ~ exponential(1);

  kappa_tau_mean ~ normal(0.985,0.1);
  kappa_omega_mean ~ normal(1.083,0.1);
  kappa_tau_sd ~ exponential(100);
  kappa_omega_sd ~ exponential(100);

  c_tau_mean ~ normal(0.001,0.001); // Informed from SS
  c_omega_mean ~ normal(0.001,0.001);
  c_tau_sd ~ exponential(50);
  c_omega_sd ~ exponential(50);

  alpha_tau_mean ~ normal(2,1);
  alpha_omega_mean ~ normal(2,1);
  alpha_tau_sd ~ exponential(30);
  alpha_omega_sd ~ exponential(50);
  
  nu_tau ~ beta(1,4);
  nu_omega ~ beta(1,4);
  
  beta_tau_omega ~ normal(1,0.05);
  beta_omega_tau ~ normal(1,0.05);
  
  rho_omega ~ normal(185,10);
  
  
  
  // Define hierarchical priors for the parameters sampled from the hierarchical distributions
  c_tau ~ normal(c_tau_mean,c_tau_sd);
  c_omega ~ normal(c_omega_mean,c_omega_sd);


  kappa_tau ~ normal(kappa_tau_mean,kappa_tau_sd);
  kappa_omega ~ normal(kappa_omega_mean,kappa_omega_sd);
  
  alpha_tau ~ normal(alpha_tau_mean,alpha_tau_sd);
  alpha_omega ~ normal(alpha_omega_mean,alpha_omega_sd);

  delta ~ normal(delta_mean,delta_sd);
  epsilon ~ normal(epsilon_mean,epsilon_sd);
  

  for (well in 1:n_wells){
    array[2*T] vector[2] mu = ode_rk45(ode_function, [c_tau[well],c_omega[well]]' .* y0[well,], 0, ts,
                                     alpha_tau[well],
                                     alpha_omega[well],
                                     // nu_tau,
                                     nu_omega,
                                     beta_tau_omega,
                                     beta_omega_tau,
                                     kappa_tau[well],
                                     kappa_omega[well]);
// need to account for biotime and geminitime
    vector[T] mu_g;
    vector[T] mu_b;
    
// Account for varying rhos.
    for (t in 1:t_star){
      mu_g[t] = rho_1*mu[t_index_g[t]][1] + rho_omega*mu[t_index_g[t]][2] + delta[well];
      mu_b[t] = mu[t_index_b[t]][1] + mu[t_index_b[t]][2] + epsilon[well];
      if (is_nan(mu_g[t])){
        mu_g[t] = 0;
      }
      if (is_nan(mu_b[t])){
        mu_b[t] = 0;
      }
    }
    for (t in (t_star+1):T){
      mu_g[t] = rho_1*mu[t_index_g[t_star]][1] + rho_omega*mu[t_index_g[t]][2] + delta[well] + rho_2*(mu[t_index_g[t]][1]-mu[t_index_g[t_star]][1]); // extra fluouresence for non tagged
      mu_b[t] = mu[t_index_b[t]][1] + mu[t_index_b[t]][2] + epsilon[well];
      if (is_nan(mu_g[t])){
        mu_g[t] = 0;
      }
      if (is_nan(mu_b[t])){
        mu_b[t] = 0;
      }
    }

    for (t in 1:T){
      y[t,well][1] ~ normal(mu_g[t],sigma * rho_2); // We scale fluourescence by rho_2
      y[t,well][2] ~ normal(mu_b[t],sigma);
    }
  }

}

