library(sumR)
library(bbmle)
source("Erlang_aux.r")
#############################
nll_adaptive <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood(data,
                          pars = theta,
                          eps =.Machine$double.eps,
                          verbose = FALSE)
}

nll_fixed <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood(data,
                          adaptive = FALSE,
                          pars = theta,
                          eps =.Machine$double.eps,
                          N_fix = 1000,
                          verbose = FALSE)
}

get_estimates <- function(data, adaptive = TRUE){
  
  inits <- log(
    c(runif(1, 1, 100),
      runif(1, .1, 50))
  )
  
  if(adaptive){
    running.time <- system.time(
      MLE <- bbmle::mle2(minuslogl = nll_adaptive,
                         start = list(log_mu = inits[1],
                                      log_beta = inits[2]),
                         method = "L-BFGS-B"
      )
    )
  }else{
    running.time <- system.time(
      MLE <- bbmle::mle2(minuslogl = nll_fixed,
                         start = list(log_mu = inits[1],
                                      log_beta = inits[2]),
                         method = "L-BFGS-B"
      )
    )
  }
  
  summy <- summary(MLE)
  
  point.ests <- as.numeric(summy@coef[, 1])
  MLE.se <- as.numeric(summy@coef[, 2])
  conf.level <- 0.95
  Z <- qnorm(p = (1 + conf.level)/2)
  
  MLE.lower <- point.ests - Z*MLE.se
  MLE.upper <- point.ests + Z*MLE.se
  
  MLE.ests <- tibble::tibble(
    parameter = c("mu", "beta"),
    true = c(Mu, B),
    point =  exp(point.ests),
    lwr = exp(MLE.lower),
    upr = exp(MLE.upper),
    time = running.time[3]
  )
  return(
    list(
      obj = MLE,
      results = MLE.ests
    )
  )
}
#############################
Mu <- 1500
B <- 1
J <- 100
data <- simulate_obsdata(n = J, mu = Mu, b = B,
                         seed = 666)

x <- data$obs_x[1]
pars <- c(Mu, B)
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(Mu, B),
                       eps = .Machine$double.eps,
                       verbose = TRUE)

est.adaptive <- get_estimates(data = data, adaptive = TRUE)

est.fixed <- get_estimates(data = data, adaptive = FALSE)

est.fixed$results
est.adaptive$results

# CIs <- exp(confint(MLE))
# CIs
# MLE.ests
## Likelihood at 'true' values:
# marginal_loglikelihood(data = data, pars = c(Mu, B),
#                            eps = .Machine$double.eps,
#                            verbose = TRUE)
## Likelihood at point estimates :
# marginal_loglikelihood(data = data, pars = c(MLE.ests$point),
#                            eps = .Machine$double.eps,
#                            verbose = TRUE) 

