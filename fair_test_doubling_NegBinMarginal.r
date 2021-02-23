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
Eps <- 1E-16
mit <- 1E5
(TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE) )

d20 <-  approx_doubling(negativeBinomial_marginalised, theta = Theta, N_start = 20,
                        n0 = obsX, epsilon = Eps, max_iter = mit)
d100 <- approx_doubling(negativeBinomial_marginalised, theta = Theta, N_start = 100,
                        n0 = obsX, epsilon = Eps, max_iter = mit)
d500 <- approx_doubling(negativeBinomial_marginalised, theta = Theta, N_start = 500,
                      n0 = obsX, epsilon = Eps, max_iter = mit)

d20$iter
d100$iter
d500$iter

robust_difference(d20$sum, TrueValue)
robust_difference(d100$sum, TrueValue)
robust_difference(d500$sum, TrueValue)

all.equal(d20$sum, TrueValue, tolerance = 0)
all.equal(d100$sum, TrueValue, tolerance = 0)
all.equal(d500$sum, TrueValue, tolerance = 0)