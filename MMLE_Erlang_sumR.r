library(sumR)
library(bbmle)
source("Erlang_aux.r")
source("Erlang_hessian.R")
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
nll_adaptive_full <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full(data,
                               pars = theta,
                               eps =.Machine$double.eps,
                               verbose = FALSE)
}
nll_fixed_full <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full(data,
                               adaptive = FALSE,
                               pars = theta,
                               eps =.Machine$double.eps,
                               N_fix = 1000,
                               verbose = FALSE)
}
nll_adaptive_full_fast <- function(log_mu, log_beta){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full_fast(data,
                                 pars = theta,
                                 eps =.Machine$double.eps,
                                 verbose = FALSE)
}

get_estimates <- function(data, minusloglik_fun,
                          conf.level = 0.95,
                          compute_Hessian = TRUE){
  
  inits <- log(
    c(runif(1, 1, 100),
      runif(1, .1, 50))
  )
  
  running.time <- system.time(
    MLE <- bbmle::mle2(minuslogl = minusloglik_fun,
                       start = list(log_mu = inits[1],
                                    log_beta = inits[2]),
                       method = "L-BFGS-B"
    )
  )
  
  summy <- summary(MLE)
  
  point.ests <- as.numeric(summy@coef[, 1])
  MLE.se <- as.numeric(summy@coef[, 2])
 
  Z <- qnorm(p = (1 + conf.level)/2)
  
  MLE.lower <- point.ests - Z*MLE.se
  MLE.upper <- point.ests + Z*MLE.se
  
  if(compute_Hessian){
    Hess <- erlangHessian(x = data$obs_x,
                          mu = exp(point.ests)[1],
                          beta = exp(point.ests)[2])
    exactVcov <- solve(-Hess)
    exactSEs <- sqrt(diag(exactVcov))
    
    exact.MLE.lower <- exp(point.ests) - Z*exactSEs
    exact.MLE.upper <- exp(point.ests) + Z*exactSEs
  }else{
    exact.MLE.lower <- rep(NA, 2)
    exact.MLE.upper <- rep(NA, 2)
  }
  
  MLE.ests <- tibble::tibble(
    parameter = c("mu", "beta"),
    true = c(Mu, B),
    point =  exp(point.ests),
    lwr = exp(MLE.lower),
    upr = exp(MLE.upper),
    lwr_exact = exact.MLE.lower,
    upr_exact = exact.MLE.upper,
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
B <- .1
J <- 10
data <- simulate_obsdata(n = J, mu = Mu, b = B,
                         seed = 666)
pars <- c(Mu, B)
##############

## Likelihood at 'true' values:
marginal_loglikelihood(data = data, pars = c(Mu, B),
                       eps = .Machine$double.eps,
                       verbose = TRUE)

marginal_loglikelihood_full(data = data, pars = c(Mu, B),
                            eps = .Machine$double.eps,
                            verbose = TRUE)

marginal_loglikelihood_full_fast(data = data, pars = c(Mu, B),
                              eps = .Machine$double.eps,
                              verbose = TRUE)

bench::mark(
  bessel_only = marginal_loglikelihood(data = data, pars = c(Mu, B),
                                       eps = .Machine$double.eps,
                                       verbose = FALSE),
  full = marginal_loglikelihood_full(data = data, pars = c(Mu, B),
                                     eps = .Machine$double.eps,
                                     verbose = FALSE),
  full_C = marginal_loglikelihood_full_fast(data = data, pars = c(Mu, B),
                                         eps = .Machine$double.eps,
                                         verbose = FALSE)
)

est.adaptive <- get_estimates(data = data, nll_adaptive)
est.fixed <- get_estimates(data = data, nll_fixed)
est.adaptive.full <- get_estimates(data = data, nll_adaptive_full)
est.fixed.full <- get_estimates(data = data, nll_fixed_full)
est.adaptive.full.fast <- get_estimates(data = data, nll_adaptive_full_fast)

est.fixed$results
est.fixed.full$results
est.adaptive$results
est.adaptive.full$results
est.adaptive.full.fast$results



##### BOOTSTRAP STUFF; gives vert narrow CIs
# ests_adaptive_bessel <- function(dt, indices){
#   dd <- dt
#   dd$obs_x <- dt$obs_x[indices]
#   obj <- get_estimates(data = dd, adaptive = TRUE)
#   return(obj$results$point)
# }
# ests_fixed_bessel <- function(dt, indices){
#   dd <- dt
#   dd$obs_x <- dt$obs_x[indices]
#   obj <- get_estimates(data = dd, adaptive = FALSE)
#   return(obj$results$point)
# }
# ests_fixed_full <- function(dt, indices){
#   dd <- dt
#   dd$obs_x <- dt$obs_x[indices]
#   obj <- get_estimates_full(data = dd, adaptive = FALSE)
#   return(obj$results$point)
# }
# ests_adaptive_full_fast <- function(dt, indices){
#   dd <- dt
#   dd$obs_x <- dt$obs_x[indices]
#   cat("data:", dt$obs_x, dd$obs_x)
#   obj <- get_estimates(data = dd, nll_adaptive_full)
#   return(obj$results$point)
# }
# 
# Nboot <- 100
# library(boot)
# system.time(
#   {
#     boot.fixed.bessel <- boot(data = data,
#                               statistic = ests_fixed_bessel,
#                               R = Nboot,
#                               parallel = "multicore",
#                               ncpus = 8)
#     boot.adaptive.bessel <- boot(data = data,
#                                  statistic = ests_adaptive_bessel,
#                                  R = Nboot,
#                                  parallel = "multicore",
#                                  ncpus = 8)
#     boot.fixed.full <- boot(data = data,
#                             statistic = ests_fixed_full,
#                             R = Nboot,
#                             parallel = "multicore",
#                             ncpus = 8)
#     boot.adaptive.full <- boot(data = data,
#                                statistic = ests_adaptive_full,
#                                R = Nboot,
#                                parallel = "multicore",
#                                ncpus = 8)
#   }
# )
# 
# fname <- paste("bootstrap_mu=", Mu,
#                "_beta=", B,
#                "_K=", J,
#                "_B=", Nboot,
#                ".RData",
#                sep = "")
# save(data,
#      boot.fixed.bessel,
#      boot.adaptive.bessel,
#      boot.fixed.full,
#      boot.adaptive.full,
#      file = fname)