functions{
  #include the_pmfs.stan
  #include infiniteSumToThreshold_1.stan
  #include infiniteSumToThreshold_2.stan
  #include compute_norm_const_adaptive.stan
  #include compute_marg_prob_adaptive.stan
}
data{
  int<lower=0> K;
  vector[K] n;
  int<lower=0> x_int[K];
  real<lower=0> x[K];
  real<lower=0> s_mu;
  real<lower=0> r_mu;
  real<lower=0> nu_sd;
  real<lower=0> a_p;
  real<lower=0> b_p;
  real<lower=0> eps;
  int<lower=0> M;
  int<lower=0, upper = 1> bayesian;
}
parameters{
  real<lower=0> mu;
  real<lower=0> phi;
  real<lower=0, upper=1> p;
}
transformed parameters{
  real log_norm_const[2] = log_Z_DOP_adaptive(mu, phi, eps, M);
  real log_marg_probs[2, K];
  for(k in 1:K){
    log_marg_probs[, k] = margProb_DOP_adaptive(mu, phi, p,
                                              x[k], x_int[k], eps, M);
  }
}
model{
  if(bayesian == 1){ // add priors
  mu ~ gamma(s_mu, r_mu);
  phi ~ normal(0, nu_sd);
  p ~ beta(a_p, b_p);
  }
  // Likelihood
  target += dot_product(n,
  to_vector(log_marg_probs[1, ]) - rep_vector(log_norm_const[1], K));
}
