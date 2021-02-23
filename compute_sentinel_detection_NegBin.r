source("truncation_package.r")
source("truncation_aux.r")

negativeBinomial_sentinel_detection <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  return(
    dnbinom(k, mu = mu, size = phi, log = TRUE) + k * log1p(-eta)
  )
}

Mu <- 20
Phi <- 1.5
Eta <- .08
Theta <- c(Mu, Phi, Eta)
Eps <- 1E-16
mit <- 1E5
(TrueValue <- Phi * (log(Phi) - log(Eta*Mu + Phi)) )

lps <- negativeBinomial_sentinel_detection(k = 0:1E4, theta = Theta)
guesslogL <- tail(diff(lps), 2)[2] 
truelogL <- log(Mu)- log_sum_exp(c(log(Mu), log(Phi))) + log1p(-Eta)

robust_difference(guesslogL, truelogL)

approx_adaptive(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = 1E5, logL = truelogL)

result <- compare_approximations(negativeBinomial_sentinel_detection, theta = Theta,
                       exact = TrueValue, epsilon = Eps, max_iter = mit, logL = truelogL)

result
# Let's see what happens when we plug in an estimate of L
compare_approximations(negativeBinomial_sentinel_detection, theta = Theta,
                       exact = TrueValue, epsilon = Eps, max_iter = 1E5, logL = guesslogL)

adapt <-  approx_adaptive(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit, logL = guesslogL)

adapt
robust_difference(TrueValue, log_sum_exp(c(adapt$Sb, adapt$St))) 
robust_difference(TrueValue, log_sum_exp(c(adapt$Sb, adapt$St.est))) ## using the estimate for the tail helps

# library(microbenchmark)
# 
# timing <- microbenchmark(
#   naive = approx_naive(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
#   naive_threshold = approx_naive_under(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
#   doubling = approx_doubling(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit),
#   adaptive = approx_adaptive(negativeBinomial_sentinel_detection, theta = Theta, epsilon = Eps, max_iter = mit)
# )
# 
# result
# timing
