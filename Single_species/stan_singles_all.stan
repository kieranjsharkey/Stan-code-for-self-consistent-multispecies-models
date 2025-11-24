functions {
  // logistic model
  vector ode_function_logistic(real t,
  vector y,
  real alpha,
  real kappa
  ) {
    vector[1] dydt;
    dydt[1] = alpha * (y[1]) * (1- (y[1]/kappa));
    return(dydt);
  }
  // Gompertz model
  vector ode_function_gompertz(real t,
  vector y,
  real alpha,
  real kappa
  ) {
    vector[1] dydt;
    dydt[1] = alpha * (y[1]) * log(kappa/y[1]);
    return(dydt);
  }
  // Richards model
  vector ode_function_richards(real t,
  vector y,
  real alpha,
  real kappa,
  real nu
  ) {
    vector[1] dydt;
    dydt[1] = alpha * (y[1]) * (1- (y[1]/kappa)^nu);
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
  int strain; // Untagged = 1, Tagged =2
  int model_type; // logistic = 1, Gompertz = 2, Richards = 3
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma;
  real<lower=0.088,upper = 0.11> epsilon_mean;
  real<lower=0> epsilon_sd;
  
  real<lower=0> c_mean;
  real<lower=0> c_sd;
  
  real<lower=0> kappa_mean;
  real<lower=0> kappa_sd;
  
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
  
  real<lower=0> nu; // If Richards, create nu.
  
  
  vector<lower=0>[n_wells] c;
  
  vector<lower=0>[n_wells] kappa;
  
  vector<lower=0>[n_wells] alpha;
  
  vector<lower=0>[n_wells] epsilon;
}
  
  // The model to be estimated. We model the output
  // 'y' to be normally distributed with mean 'mu'
  // and standard deviation 'sigma'.
  model {
    
    // Define priors
    sigma ~ exponential(100);
    epsilon_mean ~ normal(0.09555,0.001723);
    epsilon_sd ~ exponential(1);
    
    kappa_mean ~ normal(1,1);
    kappa_sd ~ exponential(1);
    
    c_mean ~ normal(0.01,0.01);
    c_sd ~ exponential(50);
    
    alpha_mean ~ normal(1,1);
    alpha_sd ~ exponential(1);
    
    if(model_type == 3) nu ~ normal(1,1);
    
    
    // Define hierarchical priors for the parameters sampled from the hierarchical distributions
    c ~ normal(c_mean,c_sd);
    
    kappa ~ normal(kappa_mean,kappa_sd);
    
    alpha ~ normal(alpha_mean,alpha_sd);

    epsilon ~ normal(epsilon_mean,epsilon_sd);

    for (well in 1:n_wells){
      array[2*T] vector[1] mu;
      if(model_type == 1){
        mu = ode_rk45(ode_function_logistic, [c[well]]', 0, ts,
        alpha[well],
        kappa[well]);
      }
      else if (model_type == 2){
        mu = ode_rk45(ode_function_gompertz, [c[well]]', 0, ts,
        alpha[well],
        kappa[well]);
      }
      else{
        mu = ode_rk45(ode_function_richards, [c[well]]', 0, ts,
        alpha[well],
        kappa[well],
        nu);
      }
      vector[T] mu_b;
      
      for (t in 1:T){
        mu_b[t] = mu[t_index_b[t]][1] + epsilon[well];
        if (is_nan(mu_b[t])){
          mu_b[t] = 0;
        }
      }
      
      
      for (t in 1:T){
        y[t,well][2] ~ normal(mu_b[t],sigma);
      }
    }
    
  }
  
  