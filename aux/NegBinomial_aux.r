negativeBinomial_marginalised_lpmf_R <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  x <- theta[4]
  ans <- ifelse(k < x,  -Inf, 
                dnbinom(k, mu = mu, size = phi, log = TRUE) +
                  dbinom(x = x, size = k, prob = eta, log = TRUE))
  return(ans)
}

Rcpp::cppFunction(code='
  NumericVector negativeBinomial_marginalised_lpmf_C(IntegerVector n, NumericVector p)
  {
  NumericVector output(n.size());
  for (int i = 0; i < n.size(); i++) {
    if (n[i] < p[3]) output[i] = -INFINITY;
    else
      output[i] = R::dnbinom_mu(n[i], p[1], p[0], true) +
        R::dbinom(p[3], n[i], p[2], true);
  }
  return output;
  }')