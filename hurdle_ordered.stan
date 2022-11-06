
data {
  
  int<lower=2> N_cutpoints;
  int<lower=0> N_obs;
  int<lower=1> P;
  array[N_obs] int y;

}

parameters {
  ordered[N_cutpoints - 1] cutpoints;
  real p;
}

model {

  for (n in 1:N_obs) {
    
    if (y[n] == 99) {
      target += bernoulli_lpmf(1 | p); 
    } else { 
      target += bernoulli_lpmf(0 | p) +  
             ordered_logistic_lpmf(y[n] | 0.0, cutpoints);
    }

  }
  
}