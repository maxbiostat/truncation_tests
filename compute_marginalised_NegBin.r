source("truncation_package.r")
source("truncation_aux.r")

negativeBinomial_marginalised <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  x <- theta[4]
  ans <- ifelse(k < x,  -Inf, 
                dnbinom(k, mu = mu, size = phi, log = TRUE) + dbinom(x = x, size = k, prob = eta, log = TRUE))
  return(ans)
}

Mu <- 50
Phi <- .5
Eta <- .008
obsX <- 20
Theta <- c(Mu, Phi, Eta, obsX)
Eps <- 1E-20
(TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE) )

result <- compare_approximations(negativeBinomial_marginalised, theta = Theta,
                       exact = TrueValue, epsilon = Eps, n0 = obsX, max_iter = 1E5)

result

library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(negativeBinomial_marginalised, theta = Theta,
                       epsilon = Eps, n0 = obsX, max_iter = mit),
  naive_threshold = approx_naive_under(negativeBinomial_marginalised, theta = Theta,
                                       epsilon = Eps, n0 = obsX, max_iter = mit),
  doubling = approx_doubling(negativeBinomial_marginalised, theta = Theta,
                             n0 = obsX, epsilon = Eps, max_iter = mit),
  adaptive = approx_adaptive(negativeBinomial_marginalised, theta = Theta,
                             n0 = obsX, epsilon = Eps, max_iter = mit)
)

result
timing