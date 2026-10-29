#include functions.stan
data {
 int<lower = 0> N;         // number of points
 vector[N] y;              // data
}

transformed data {
  real median = quantile(y, 0.5);              // Median of y
  
  vector[N] abs_dev;
  for (n in 1:N) {
    abs_dev[n] = abs(y[n] - median);
  }
  real mad = quantile(abs_dev, 0.5) * 1.4826;  // MAD scaled by 1.4826 for consistency with R
  
  vector[N] y_std = (y - median) / mad;        // Standardised y data
}

parameters {
  real<lower=0.8, upper=0.95> frac_neut;    // Fraction of density that is neutral
  //real<lower=1e-3, upper=0.2> sigma_neut; // Standard deviation of the neutral peak
  real<lower=0.9, upper=2.5> beta;         // Scale
  real<lower=-0.6> delta;                   // Location
}

transformed parameters {
  real<lower=0, upper=1> theta_noise;
  real<lower=0, upper=1> theta_ko;
  real<lower=0, upper=1> theta_del;
  real<lower=0, upper=1> theta_neut;
  real<lower=0, upper=1> theta_pos;
  real<lower=0, upper=1> sigma_neut;
  
  sigma_neut = 0.15;
  
  theta_noise = 0.03;
  theta_ko = 0.001;
  theta_del = fmax(0.5 * (1 - frac_neut - theta_noise - theta_ko), 1e-6);
  theta_neut = fmax(frac_neut, 1e-6);;
  theta_pos = fmax(0.5 * (1 - frac_neut - theta_noise - theta_ko), 1e-6);
}

model {
  frac_neut ~ beta(3, 1);
  
  // sigma_neut ~ lognormal(log(0.15), 0.04);
  
  beta ~ lognormal(log(1), 0.01);
  delta ~ normal(0, 0.15); 
  
  for (n in 1:N) {
    target += model_log_prob(f_inverse(y_std[n], beta, delta), sigma_neut, theta_ko, theta_del, theta_neut, theta_pos, theta_noise);
  }
}

generated quantities {
  real y_median = median;          // Median value of y
  real y_mad = mad;                // MAD of y
  vector[N] log_lik;               // Per data point Log-likelihood
  
  for (n in 1:N) {
    log_lik[n] = model_log_prob(f_inverse(y_std[n], beta, delta), sigma_neut, theta_ko, theta_del, theta_neut, theta_pos, theta_noise);
  }
}
