source("truncation_package.r")
source("truncation_aux.r")

weird_series <- function(k, theta){
  s <- theta[1]
  if(k==0) return(-Inf)
  return(
    -(2*log(k) + k*log(s) )
  )
}
weird_series <- Vectorize(weird_series)

S <- 2
Theta <- S
Eps <- 1E-20
TrueValue <- log(copula::polylog(z = 1/S, s = 2))
lgL <- -log(S)

result <- compare_approximations(weird_series, theta = Theta,
                                 exact = TrueValue, epsilon = Eps, max_iter = 1E5, logL = lgL)

library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  naive_threshold = approx_naive_under(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  doubling = approx_doubling(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  adaptive = approx_adaptive(weird_series, theta = Theta, epsilon = Eps, max_iter = mit)
)

result
timing