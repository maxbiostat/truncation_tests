real lbb(real x, int n, real a, real b){
  if(n < x){
    return negative_infinity();
  }else{
    return lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b);
  }
}
real logFunction(int k, real[] pars){
  real mu = pars[1];
  real a = pars[2];
  real b = pars[3];
  real x = pars[4];
  real ans =  lbb(x, k, a, b) + poisson_lpmf(k | mu);
  return(ans);
}
