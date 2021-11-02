functions {
  #include R0_example_pmfs.stan
  #include infiniteAdaptive.stan
  real compute_log_L(real[] theta){
    real r = theta[1];
    real w = theta[2];
    real ans = (1 + w) * ( log(1 + w) - log(w + r) ) + log(r) ;
    return(ans);
  }
  real[] log_marg_prob_s(int s,
                         real r, real w,
                         real a, real b,
                      real epsilon, int max_it, real logL){
    int x0;
    if(s == 0){
      x0 = 1;
    }else{
      x0 = s;
    }
    real lprob[2] = infiniteAdaptive({r, w, a , b, s}, epsilon, max_it, logL, x0);
    return(lprob);
  }
}
data {
  int<lower=1> D; // number of unique sequence cluster sizes
  int<lower=1> s[D]; // sequence cluster sizes (censored)
  vector[D] n; // size counts
  real<lower=0> epsilon;
  int<lower=0> max_iter;
}
parameters {
  real<lower=0, upper=1> R0;
  real<lower=0> omega;
  real<lower=0> alpha;
  real<lower=0> beta;
}
transformed parameters{
  real log_Prs[2, D];
  real lgL = compute_log_L({R0, omega});
  real lp0[2] = log_marg_prob_s(0, R0, omega, alpha, beta,
  epsilon, max_iter, lgL);
  for(i in 1:D){
    log_Prs[, i] =  log_marg_prob_s(s[i], R0, omega,
                                     alpha, beta, 
                                     epsilon, max_iter, lgL);
  } 
}
model{
   target += dot_product(n,
                      to_vector(log_Prs[1, ]) - rep_vector(log1m_exp(lp0[1]), D));
}
generated quantities{
  real<lower=0, upper=1> p0 = exp(lp0[1]);
}
