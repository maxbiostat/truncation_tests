source("truncation_package.r")
source("truncation_aux.r")

poisson_sentinel_detection <- function(k, theta){
  lambda <- theta[1]
  eta <- theta[2]
  return(
    dpois(k, lambda = lambda, log = TRUE) + k * log(1-eta)
  )
}

Lambda <- 200
Eta <- .08
Theta <- c(Lambda, Eta)
Eps <- 1E-20
(TrueValue <- -Lambda * Eta)

result <- compare_approximations(poisson_sentinel_detection, theta = Theta,
                       exact = TrueValue, epsilon = Eps, max_iter = 1E5)


library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(poisson_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
  naive_threshold = approx_naive_under(poisson_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
  doubling = approx_doubling(poisson_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
  adaptive = approx_adaptive(poisson_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit)
)

result
timing
