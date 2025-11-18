functions {
  vector ode_function(real t,
                      vector y,
                      real alpha_A,
                      real alpha_B,
                      real beta_AB,
                      real beta_BA,
                      real kappa_A,
                      real kappa_B
  ) {
  vector[2] dydt;
  dydt[1] = alpha_A * (y[1]) * (1 - (y[1] / kappa_A) - beta_AB * (y[2] / kappa_B));
  dydt[2] = alpha_B * (y[2]) * (1 - (y[2] / kappa_B) - beta_BA * (y[1] / kappa_A));
    
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
  real sigma_A;
  real rho_1;
  real rho_2;
  int tau;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma_B;
  real<lower=0,upper=380> epsilon_A_mean; // Was 270
  real<lower=0,upper = 0.11> epsilon_B_mean;
  real<lower=0> epsilon_A_sd;
  real<lower=0> epsilon_B_sd;
  real<lower=0> B0_A_mean;
  real<lower=0> B0_B_mean;
  
  real<lower=130, upper=240> rho_est;

  // real<lower=0> B0_mean;

  real<lower=0> kappa_A_mean;
  real<lower=0> kappa_B_mean;
  real<lower=0> B0_A_sd;
  real<lower=0> B0_B_sd;

  // real<lower=0> B0_sd;

  real<lower=0> kappa_A_sd;
  real<lower=0> kappa_B_sd;
  real<lower=0> alpha_A_mean;
  real<lower=0> alpha_B_mean;
  real<lower=0> alpha_A_sd;
  real<lower=0> alpha_B_sd;

  real<lower=0> beta_AB;
  real<lower=0> beta_BA;
  vector<lower=0>[n_wells] B0_A;
  vector<lower=0>[n_wells] B0_B;

  vector<lower=0>[n_wells] kappa_A;
  vector<lower=0>[n_wells] kappa_B;
  
  vector<lower=0>[n_wells] alpha_A;
  vector<lower=0>[n_wells] alpha_B;
  
  vector<lower=0,upper=380>[n_wells] epsilon_A;
  vector<lower=0,upper = 0.11>[n_wells] epsilon_B;
}

// transformed parameters {
//   array[2] vector[n_wells] B0_vect = {B0_A, B0_B};
// }

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // Define priors
  sigma_B ~ exponential(100);
  epsilon_A_mean ~ normal(331.16,12.94);
  epsilon_B_mean ~ normal(0.09555,0.001723);
  epsilon_A_sd ~ exponential(100);
  epsilon_B_sd ~ exponential(100);

  kappa_A_mean ~ normal(0.985,0.1);
  kappa_B_mean ~ normal(1.083,0.1);
  kappa_A_sd ~ exponential(1000);
  kappa_B_sd ~ exponential(1000);

  B0_A_mean ~ normal(0.01,0.01); // Informed from SS
  B0_B_mean ~ normal(0.01,0.01);
  B0_A_sd ~ exponential(50);
  B0_B_sd ~ exponential(50);

  alpha_A_mean ~ normal(0.8,1);
  alpha_B_mean ~ normal(0.8,1);
  alpha_A_sd ~ exponential(30);
  alpha_B_sd ~ exponential(50);
  
  beta_AB ~ normal(1,0.05);
  beta_BA ~ normal(1,0.05);
  
  rho_est ~ normal(185,10);
  
  // Define hierarchical priors for the parameters sampled from the hierarchical distributions
  B0_A ~ normal(B0_A_mean,B0_A_sd);
  B0_B ~ normal(B0_B_mean,B0_B_sd);

  kappa_A ~ normal(kappa_A_mean,kappa_A_sd);
  kappa_B ~ normal(kappa_B_mean,kappa_B_sd);
  
  alpha_A ~ normal(alpha_A_mean,alpha_A_sd);
  alpha_B ~ normal(alpha_B_mean,alpha_B_sd);

  epsilon_A ~ normal(epsilon_A_mean,epsilon_A_sd);
  epsilon_B ~ normal(epsilon_B_mean,epsilon_B_sd);
  

  // Want to force initial distirbution to be out of 100%?
  for (well in 1:n_wells){
    array[2*T] vector[2] mu = ode_rk45(ode_function, [B0_A[well],B0_B[well]]' .* y0[well,], 0, ts,
                                     alpha_A[well],
                                     alpha_B[well],
                                     beta_AB,
                                     beta_BA,
                                     kappa_A[well],
                                     kappa_B[well]);
                  // need to account for biotime and geminitime
    vector[T] mu_g;
    vector[T] mu_b;
// Edit epsilons.
    for (t in 1:tau){
      mu_g[t] = rho_1*mu[t_index_g[t]][1] + rho_est*mu[t_index_g[t]][2] + epsilon_A[well];
      mu_b[t] = mu[t_index_b[t]][1] + mu[t_index_b[t]][2] + epsilon_B[well];
      if (is_nan(mu_g[t])){
        mu_g[t] = 0;
      }
      if (is_nan(mu_b[t])){
        mu_b[t] = 0;
      }
    }
    for (t in (tau+1):T){
      mu_g[t] = rho_1*mu[t_index_g[tau]][1] + rho_est*mu[t_index_g[t]][2] + epsilon_A[well] + rho_2*(mu[t_index_g[t]][1]-mu[t_index_g[tau]][1]); // extra fluouresence for non tagged
      mu_b[t] = mu[t_index_b[t]][1] + mu[t_index_b[t]][2] + epsilon_B[well];
      if (is_nan(mu_g[t])){
        mu_g[t] = 0;
      }
      if (is_nan(mu_b[t])){
        mu_b[t] = 0;
      }
    }

    for (t in 1:T){
      y[t,well][1] ~ normal(mu_g[t],sigma_B * sigma_A); // test with a different variance for fluoresence
      y[t,well][2] ~ normal(mu_b[t],sigma_B);
    }
  }

}

