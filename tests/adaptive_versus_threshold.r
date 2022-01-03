# remotes::install_github("GuidoAMoreira/sumR")
library(sumR)

Rcpp::cppFunction(code='
  double compute_lterm(long n, Rcpp::NumericVector p)
  {
  long double ans = - 2*log(n+1) - (n+1)* log(p[0]);
  return ans;
  }')

my_dilog <- function(z) 
{
  if (abs(z) >= 1) {
    return(NA)
  }
  out = stats::integrate(function(y) log1p(-y)/y, lower = 0, 
                         upper = z, subdivisions = 1E3,
                         rel.tol = 1E4 * .Machine$double.eps)$value
  return(-out)
}

##############3

Eps <- .Machine$double.eps *10
M <- 1E5
a1 <- 2

lps.1 <- sapply(0:(2*M),
                function(j) compute_lterm(n = j, p = a1))
TV.1 <- matrixStats::logSumExp(sort(lps.1))
TTV.1 <- log(1/12 * (pi^2 - 6*log(2)^2))

naive.1 <- sumR::infiniteSum(
  logFunction = compute_lterm,
  parameters = a1,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  forceAlgorithm = 1
)

adaptive.1 <- sumR::infiniteSum(
  logFunction = compute_lterm,
  parameters = a1,
  epsilon = Eps,
  maxIter = M,
  logL = -log(a1),
  n0 = 0,
  forceAlgorithm = 2
)

doubling.1 <- sumR::infiniteSum_cFolding_C(
  logFunction = compute_lterm,
  parameters = a1,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  c = 2,
  N_start = 10
)

####################
a2 <- 1.0005
L2 <- 1/a2

lps.2 <- sapply(0:(2*M),
                function(j) compute_lterm(n = j, p = a2))
TV.2 <- matrixStats::logSumExp(sort(lps.2))
TTV.2 <- log(my_dilog(1/a2))

naive.2 <- sumR::infiniteSum(
  logFunction = compute_lterm,
  parameters = a2,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  forceAlgorithm = 1
)

adaptive.2 <- sumR::infiniteSum(
  logFunction = compute_lterm,
  parameters = a2,
  epsilon = Eps,
  maxIter = M,
  logL = -log(a2),
  n0 = 0,
  forceAlgorithm = 2
)

doubling.2 <- sumR::infiniteSum_cFolding_C(
  logFunction = compute_lterm,
  parameters = a2,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  c = 2,
  N_start = 10
)

source("../aux/aux.r")

robust_difference(x = TV.1, y = TTV.1)

robust_difference(x = TV.1, y = naive.1$sum)
robust_difference(x = TV.1, y = adaptive.1$sum)
robust_difference(x = TV.1, y = doubling.1$sum)
Eps

robust_difference(x = TTV.1, y = naive.1$sum)
robust_difference(x = TTV.1, y = adaptive.1$sum)
robust_difference(x = TTV.1, y = doubling.1$sum)
Eps

e2 <- Eps * exp(-log(a2)-log1p(-L2))  ## This should be theoretical bound

robust_difference(x = TV.2, y = TTV.2)

robust_difference(x = TV.2, y = naive.2$sum)
robust_difference(x = TV.2, y = adaptive.2$sum)
robust_difference(x = TV.2, y = doubling.2$sum)
Eps
e2

robust_difference(x = TTV.2, y = naive.2$sum)
robust_difference(x = TTV.2, y = adaptive.2$sum)
robust_difference(x = TTV.2, y = doubling.2$sum)
Eps
e2

naive.1$n
adaptive.1$n
doubling.1$n

naive.2$n
adaptive.2$n
doubling.2$n