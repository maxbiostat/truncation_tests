map_Rn <- function(x){
  return(
    c(
      arm::logit(x[1]),
      log(x[2]),
      log(x[3]),
      log(x[4])
    )
  )
}
map_back <- function(y){
  return(
    c(
      arm::invlogit(y[1]),
      exp(y[2]),
      exp(y[3]),
      exp(y[4])
    )
  )
}
##
library(sumR)
library(extraDistr)
library(stats4)
source("R0_Betabinomial_truc_aux.r")
####################
R0 <- .67
Disp <- .5
A <- 1
B <- 1
J <- 1000
simu <- simulate_obsdata(n = J,
                         r = R0, w = Disp,
                         a = A, b = B,
                         seed = 666)
data <- compress_counts(simu$obs_xprime)

mean(simu$obs_y)
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(R0, Disp, A, B),
                       verbose = TRUE)

####

nll_adaptive <- function(logit_R0, log_omega, log_alpha, log_beta){
  theta <- map_back(c(logit_R0, log_alpha, log_beta))
  -marginal_loglikelihood(data, pars = theta)
}

inits <- map_Rn(c(.5, .5, 1, 1))

MLE <- bbmle::mle2(minuslogl = nll_adaptive,
                       start = list(logit_R0 = inits[1],
                                    log_omega = inits[2],
                                    log_alpha = inits[3],
                                    log_beta = inits[4]),
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
  true = c(R0, Disp, A, B),
  mean =  map_back(point.ests),
  lwr = map_back(MLE.lower),
  upr = map_back(MLE.upper)
)
MLE.ests
# CIs <- confint(MLE)
