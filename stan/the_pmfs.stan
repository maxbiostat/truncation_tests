/* Double Poisson probability mass function */ 
real log_DOP(int k, real mu, real phi){
  if(k==0){
    return(negative_infinity());
  }else{
    return(0.5 * log(phi) - (mu - k) * phi - k +
    k * (1 - phi) * log(k) - lgamma(k + 1) + phi * k * log(mu));
  }
}
real logFunction(int n, real[] p){
  return(log_DOP(n, p[1], p[2]));
}

/* Marginalised (unnormalised) probability mass function */ 

real bin_log(real x, real n, real p){
  if(n < x){
    return negative_infinity();
  }else{
    return lchoose(n, x) + x*log(p) + (n-x)*log1m(p);  
  }
}
real log_marg(int k, real mu, real phi, real p, real x){
  return log_DOP(k, mu, phi) + bin_log(x, k, p);
}
real logFunction_2(int n, real[] p){
  return(log_marg(n, p[1], p[2], p[3], p[4]));
}
