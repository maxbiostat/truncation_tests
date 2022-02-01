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

Eps <- .Machine$double.eps
M <- 5E5
a1 <- 2
L1 <- 1/a1
B1 <- 1/(a1-1)

lps.1 <- sapply(0:(1E6),
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

batches.1 <- sumR::infiniteSum_batches_C(
  logFunction = compute_lterm,
  parameters = a1,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  batch_size = 2
)

####################
a2 <- 1.1
L2 <- 1/a2
B2 <- 1/(a2-1)

lps.2 <- sapply(0:(1E6),
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

batches.2 <- sumR::infiniteSum_batches_C(
  logFunction = compute_lterm,
  parameters = a2,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  batch_size = round(B2/2)
)

batches.2.correct <- sumR::infiniteSum_batches_C(
  logFunction = compute_lterm,
  parameters = a2,
  epsilon = Eps,
  maxIter = M,
  n0 = 0,
  batch_size = B2 + 1
)

source("../aux/aux.r")

out1a <- tibble::tibble(
  method = c("Threshold", "BoundingPair", "Batches"),
  target_error = Eps,
  L = L1,
  a = a1,
  error = c(
    robust_difference(x = TV.1, y = naive.1$sum),
    robust_difference(x = TV.1, y = adaptive.1$sum),
    robust_difference(x = TV.1, y = batches.1$sum)
  ),
  n_iter = c(
    naive.1$n,
    adaptive.1$n,
    batches.1$n
  ),
  true_method = "Huge"
)

out1b <- tibble::tibble(
  method = c("Threshold", "BoundingPair", "Batches"),
  target_error = Eps,
  L = L1,
  a = a1,
  error = c(
    robust_difference(x = TTV.1, y = naive.1$sum),
    robust_difference(x = TTV.1, y = adaptive.1$sum),
    robust_difference(x = TTV.1, y = batches.1$sum)
  ),
  n_iter = c(
    naive.1$n,
    adaptive.1$n,
    batches.1$n
  ),
  true_method = "Closed-form"
)

out1 <- rbind(out1a, out1b)


out2a <- tibble::tibble(
  method = c("Threshold", "BoundingPair",
             "Batches_wrong", "Batches_right"),
  target_error = Eps,
  L = L2,
  a = a2,
  error = c(
    robust_difference(x = TV.2, y = naive.2$sum),
    robust_difference(x = TV.2, y = adaptive.2$sum),
    robust_difference(x = TV.2, y = batches.2$sum),
    robust_difference(x = TV.2, y = batches.2.correct$sum)
  ),
  n_iter = c(
    naive.2$n,
    adaptive.2$n,
    batches.2$n,
    batches.2.correct$n
  ),
  true_method = "Huge"
)

out2b <- tibble::tibble(
  method = c("Threshold", "BoundingPair",
             "Batches_wrong", "Batches_right"),
  target_error = Eps,
  L = L2,
  a = a2,
  error = c(
    robust_difference(x = TTV.2, y = naive.2$sum),
    robust_difference(x = TTV.2, y = adaptive.2$sum),
    robust_difference(x = TTV.2, y = batches.2$sum),
    robust_difference(x = TTV.2, y = batches.2.correct$sum)
  ),
  n_iter = c(
    naive.2$n,
    adaptive.2$n,
    batches.2$n,
    batches.2.correct$n
  ),
  true_method = "dilogarithm"
)

out2 <- rbind(out2a, out2b)

##################

out1 
robust_difference(x = TV.1, y = TTV.1)

out2
Eps * exp(-log(a2)-log1p(-L2))  ## This should be theoretical upper bound on the error
robust_difference(x = TV.2, y = TTV.2)

subset(rbind(out1, out2), true_method == "Huge")
B1
B2


