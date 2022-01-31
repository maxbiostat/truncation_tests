poisson_lfactmom_R <- function(k, theta){
  x <- theta[2]
  ans <- ifelse(k < x, -Inf,
                lfactorial(k) - lfactorial(k-x) +
                  dpois(x = k, lambda = theta[1], log = TRUE))
  return(ans)
}

Rcpp::cppFunction(code='
  NumericVector poisson_lfactmom_C(IntegerVector n, NumericVector p)
  {
  NumericVector output(n.size());
  for (int i = 0; i < n.size(); i++) {
    if (p[1] > n[i]) output[i] = -INFINITY;
    else
      output[i] = R::dpois(n[i], p[0], true) +
      lgamma(n[i] + 1) - lgamma(n[i] - p[1] + 1);
  }
  return output;
  }')