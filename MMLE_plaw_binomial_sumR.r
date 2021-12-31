library(sumR)
library(extraDistr)
library(stats4)
source("plaw_binomial_trunc_aux.r")

map_Rn <- function(x){
  return(
    c(
      log(x[1]),
      arm::logit(x[2])
    )
  )
}
map_back <- function(y){
  return(
    c(
      exp(y[1]),
      arm::invlogit(y[2])
    )
  )
}

####################
Xm <- 1
Alpha <- 1.5
detProb <- .01
J <- 500
simu <- simulate_obsdata(n = J,
                         x_m = Xm,
                         alpha = Alpha,
                         p = detProb,
                         seed = 666)

data <- compress_counts(simu$obs_xprime)

mean(simu$obs_xprime)
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(Xm, Alpha, detProb),
                       verbose = TRUE)

####

lik_adaptive <- function(log_alpha, logit_p){
  theta <- c(Xm, map_back(c(log_alpha, logit_p)))
  -marginal_loglikelihood(data, pars = theta)
}

inits <- map_Rn(c(runif(1, 1.1, 10), runif(1, .01, .99)))

MLE <- bbmle::mle2(minuslogl = lik_adaptive,
                       start = list(log_alpha = inits[1],
                                    logit_p = inits[2]),
                       lower = c(log_alpha = log(1.1),
                                 logit_p = arm::logit(.01)),
                       upper = c(log_alpha = log(10),
                                 logit_p = arm::logit(.99)),
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
  parameter = c("alpha", "beta"),
  true = c(Alpha, detProb),
  point =  map_back(point.ests),
  lwr = map_back(MLE.lower),
  upr = map_back(MLE.upper)
)
MLE.ests
# CIs <- confint(MLE)
