functions {
  real f_inverse(real s, real beta, real delta) {
    // Transform applied to s to get w - could also be considered the inverse 
    // of the function s = f(w) which is more "forward" direction
    // w = 2 / (1 + exp(-(beta * (s + delta))))
    // Avoid under/overflow in exp
    real arg = fmin(fmax(beta * (s + delta), -20), 20);
    return 2 / (1 + exp(-arg));
  }
  
  real model_log_prob(real w, real sigma_neut, real theta_ko, real theta_del, real theta_neut, real theta_pos, real theta_noise) {
    vector[5] log_probs;
    
    log_probs[1] = log(theta_ko) + normal_lpdf(w | 0.1, 0.05);              // KO variants
    log_probs[2] = log(theta_del) + skew_normal_lpdf(w | 0.75, 0.7, -20);   // Deleterious variants
    log_probs[3] = log(theta_neut) + normal_lpdf(w | 1, sigma_neut);        // Neutral variants
    log_probs[4] = log(theta_pos) + skew_normal_lpdf(w | 1.25, 0.7, 20);    // Positive variants
    log_probs[5] = log(theta_noise) + beta_lpdf(w / 2 | 1, 1) - log(2);     // Noise soaking
    
    return log_sum_exp(log_probs);
  }
}