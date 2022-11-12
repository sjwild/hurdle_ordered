
data {
  
  int<lower=2> N_levels;
  int<lower=0> N_obs;
  int<lower=1> P;
  array[N_obs] int y;
  matrix[N_obs, P] X;
  int DK; 

}

parameters {
  ordered[N_levels - 1] cutpoints;
  vector[P] beta;
  vector[P] gamma;
  //real<lower=0, upper=1> p;
}

transformed parameters {

  
}

model {

  // priors
  //p ~ beta(2, 2);
  beta ~ normal(0, 1.5);
  gamma ~ normal(0, 1.5);
  
  // model
  {
  vector[N_obs] eta;
  vector[N_obs] mu;
  
  eta = X * beta;
  mu = X * gamma;
  for (n in 1:N_obs) {
    
    if (y[n] == DK) {
      target += bernoulli_logit_lpmf(1 | mu[n]); 
    } else { 
      target += bernoulli_logit_lpmf(0 | mu[n]) +  
             ordered_logistic_lpmf(y[n] | eta[n], cutpoints);
    }

  }
  
  }
  
}