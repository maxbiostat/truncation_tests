library(sumR)
library(numDeriv)
source("aux/aux.r")
source("Erlang_aux.r")
source("Erlang_hessian.R")
###########
Mu <- 1500
B <- .56
J <- 5
data <- simulate_obsdata(n = J, mu = Mu, b = B,
                         seed = 1234)

#########
bessel_lik <- function(x){
  marginal_loglikelihood(data = data,
                         pars = x,
                         eps = .Machine$double.eps,
                         verbose = FALSE)
}

full_lik <- function(x){
  marginal_loglikelihood_full_fast(data = data, pars = x,
                                   eps = .Machine$double.eps,
                                   verbose = FALSE)
}

theta_0 <- c(Mu, B+B/2)

bessel_lik(theta_0)
full_lik(theta_0)

(M1 <- numDeriv::hessian(
  bessel_lik,
  x = theta_0
) )

(M2 <- numDeriv::hessian(
  full_lik,
  x = theta_0
))

(M3 <- erlangHessian(x = data$obs_x,
                     mu = theta_0[1],
                     beta = theta_0[2]) )

M3/M2

matrixcalc::is.negative.semi.definite(M1)
matrixcalc::is.negative.semi.definite(M2)
matrixcalc::is.negative.semi.definite(M3)

sqrt(solve(-M3))
