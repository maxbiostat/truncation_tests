library(sumR)
library(extraDistr)
library(stats4)
source("poisson_Betabinomial_aux.r")
####################
Mu <- 40
A <- 1
B <- 10
J <- 1000
simu <- simulate_obsdata(n = J, mu = Mu, a = A, b = B,
                         seed = 666)
data <- compress_counts(simu$obs_x)

mean(simu$obs_y)
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(Mu, A, B),
                       verbose = TRUE)

## Profile likelihood: Mu
f_mu <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(x, A, B))
}
f_mu <- Vectorize(f_mu)
mus <- seq(from = 0.1*Mu, to = 10*Mu, length.out = 20)
plot(mus, f_mu(mus))
abline(v = Mu, lwd = 2, lty = 2)

## Profile likelihood: Phi
f_alpha <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(Mu, x, B))
}
f_alpha <- Vectorize(f_alpha)
As <- seq(from = 0.1*A, to = 5*A, length.out = 20)
plot(As, f_alpha(As))
abline(v = A, lwd = 2, lty = 2)

## Profile likelihood: p
f_beta <- function(x){
  marginal_loglikelihood(data = data, 
                         pars = c(Mu, A, x))
}
f_beta <- Vectorize(f_beta)
Bs <- seq(from = 0.1*B, to = 5*B, length.out = 20)
plot(Bs, f_beta(Bs))
abline(v = B, lwd = 2, lty = 2)

####

lik_adaptive <- function(log_mu, log_alpha, log_beta){
  theta <- exp(c(log_mu, log_alpha, log_beta))
  -marginal_loglikelihood(data, pars = theta)
}

inits <- log(c(runif(1, 1, 100),
               runif(1, 1, 100),
               runif(1, 1, 100) )
)

MLE <- bbmle::mle2(minuslogl = lik_adaptive,
                       start = list(log_mu = inits[1],
                                    log_alpha = inits[2],
                                    log_beta = inits[3]),
                       method = "L-BFGS-B"
)

summy <- summary(MLE)

point.ests <- as.numeric(summy@coef[, 1])
MLE.se <- as.numeric(summy@coef[, 2])
conf.level <- 0.95
Z <- qnorm(p = (1 + conf.level)/2)

MLE.lower <- point.ests - Z*MLE.se
MLE.upper <- point.ests + Z*MLE.se

MLE.ests <- tibble::tibble(
  parameter = c("mu", "alpha", "beta"),
  true = c(Mu, A, B),
  mean =  exp(point.ests),
  lwr = exp(MLE.lower),
  upr = exp(MLE.upper)
)
MLE.ests
# CIs <- confint(MLE)
