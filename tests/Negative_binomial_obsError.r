library(sumR)
library(tidyverse)
source("../aux/aux.r")
source("../aux/NegBinomial_aux.r")

Mu <- 100
Phi <- .1
Eta <- .01
obsX <- 0
Theta <- c(Mu, Phi, Eta, obsX)
lgL <- log(Mu) - matrixStats::logSumExp(c(log(Mu), log(Phi))) + log1p(-Eta)
L <- exp(lgL)
B <- L/(1-L)
Eps <- .Machine$double.eps
M <- 5E5
TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE)
#######
negativeBinomial_marginalised_lpmf_R(k = obsX, theta = Theta)
negativeBinomial_marginalised_lpmf_C(n = obsX, p = Theta)
#######

result.R <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_R,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   max_iter = M,
                                   logL = lgL)

result.R2 <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_R,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   max_iter = M,
                                   batch_size = B+1,
                                   logL = lgL)

result.C <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   max_iter = M,
                                   logL = lgL)

result.C2 <- compare_approximations(compute_lterm = negativeBinomial_marginalised_lpmf_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   max_iter = M,
                                   batch_size = B+1,
                                   logL = lgL)
result.R
result.R2
result.C
result.C2

( error.bound <- Eps * exp(lgL - log1m_exp(lgL)) ) ## This should be theoretical bound

