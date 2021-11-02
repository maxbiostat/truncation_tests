real size_lpmf(int y, real r, real omega) {
  real a = r/omega;
  real b = a/(1 + a)^(1 + omega);
  real lp = log(1+a)-log(a) + lchoose(omega*y + y - 2, omega*y -1) -log(y) + y*log(b);
  return lp;
}
real lbb(real x, int n, real a, real b){
  if(n < x){
    return negative_infinity();
  }else{
    return lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b);
  }
}
real logFunction(int k, real[] pars){
  real r = pars[1];
  real w = pars[2];
  real a = pars[3];
  real b = pars[4];
  real s = pars[5];
  real ans =  lbb(s, k, a, b) + size_lpmf(k | r, w);
  return(ans);
}
