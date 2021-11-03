functions {
  #include PBB_pmfs.stan
  #include infiniteSumToThreshold.stan
  real[] log_marg_prob_s(int x,
                         real mu,
                         real a, real b,
                         real epsilon, int max_it){
    int x0;
    if(x == 0){
      x0 = 1;
    }else{
      x0 = x;
    }
    real lprob[2] = infiniteSumToThreshold({mu, a , b, x}, epsilon, max_it, x0);
    return(lprob);
  }
}
data {
  int<lower=1> K; // number of unique sequence cluster sizes
  int<lower=1> x[K]; // sequence cluster sizes (censored)
  vector[K] n; // size counts
  real<lower=0> epsilon;
  int<lower=0> max_iter;
}
parameters {
  real<lower=0> mu;
  real<lower=0> alpha;
  real<lower=0> beta;
}
transformed parameters{
  real log_Prs[2, K];
  real lp0[2] = log_marg_prob_s(0, mu, alpha, beta, epsilon, max_iter);
  for(i in 1:K){
    log_Prs[, i] = log_marg_prob_s(x[i], mu, alpha, beta, 
                                     epsilon, max_iter);
  } 
}
model{
   target += dot_product(n,
                      to_vector(log_Prs[1, ]) - rep_vector(log1m_exp(lp0[1]), K));
}
generated quantities{
  real<lower=0, upper=1> p0 = exp(lp0[1]);
}
