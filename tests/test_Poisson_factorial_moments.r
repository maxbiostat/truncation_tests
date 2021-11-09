library(sumR)
source("../aux/aux.r")
#######
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

#######
Mu <- 12 ## try 11.8 and 11.9 and be amazed
R <- 2
Theta <- c(Mu, R)
TrueValue <- R*log(Mu)
Eps <- .Machine$double.eps
lgL <- log(0)
#######

result.R <- compare_approximations(compute_lterm = poisson_lfactmom_R,
                                 theta = Theta,
                                 exact = TrueValue,
                                 eps = Eps,
                                 logL = lgL)

result.C <- compare_approximations(compute_lterm = poisson_lfactmom_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)
result.R
result.C