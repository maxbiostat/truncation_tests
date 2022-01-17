library(sumR)
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
Mu <- 11.8 ## try 11.8 and 11.9 and be amazed
R <- 3
theta <- c(Mu, R)
eps <- .Machine$double.eps 
logL <- log(0)
max_iter <- 3E5
n0 <- 0
Nstart <- 20



timings <- bench::mark(
  adaptive_R = sumR::infiniteSum(
    logFunction = poisson_lfactmom_C,
    parameters = theta,
    epsilon = eps,
    logL = log(0),
    maxIter = max_iter,
    n0 = n0,
    forceAlgorithm = 2
  ),
  naive_R = sumR::infiniteSum(
    logFunction = poisson_lfactmom_C,
    parameters = theta,
    epsilon = eps,
    maxIter = max_iter,
    n0 = n0,
    forceAlgorithm = 1
  ),
  doubling_R = sumR::infiniteSum_cFolding(
    logFunction = poisson_lfactmom_C,
    parameters = theta,
    epsilon = eps,
    N_start = Nstart,
    c = 2,
    maxIter = max_iter,
    n0 = n0
  ),
  check = FALSE,
  iterations = 200
)


timings[, c("expression", "min", "median",
            "itr/sec", "gc/sec")]


adaptive_R <- sumR::infiniteSum(
  logFunction = poisson_lfactmom_C,
  parameters = theta,
  epsilon = eps,
  logL = log(0),
  maxIter = max_iter,
  n0 = n0,
  forceAlgorithm = 2
)
naive_R <- sumR::infiniteSum(
  logFunction = poisson_lfactmom_C,
  parameters = theta,
  epsilon = eps,
  maxIter = max_iter,
  n0 = n0,
  forceAlgorithm = 1
)
doubling_R <- sumR::infiniteSum_cFolding(
  logFunction = poisson_lfactmom_C,
  parameters = theta,
  epsilon = eps,
  N_start = Nstart,
  c = 2,
  maxIter = max_iter,
  n0 = n0
)
