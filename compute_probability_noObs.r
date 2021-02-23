source("truncation_package.r")
source("truncation_aux.r")

probability_noObs <- function(k, theta){
  eta <- theta[1]
  return(
    k * log(1-eta)  
  )
}

Eta <- .08
Theta <- Eta
Eps <- 1E-14
(TrueValue <- -log(Eta))

result <- compare_approximations(probability_noObs, theta = Theta,
                                 exact = TrueValue, epsilon = Eps, max_iter = 1E5)

result

library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(probability_noObs, theta = Theta, epsilon = Eps, max_iter = mit),
  naive_threshold = approx_naive_under(probability_noObs, theta = Theta, epsilon = Eps, max_iter = mit),
  doubling = approx_doubling(probability_noObs, theta = Theta, epsilon = Eps, max_iter = mit),
  adaptive = approx_adaptive(probability_noObs, theta = Theta, epsilon = Eps, max_iter = mit)
)

result
timing
