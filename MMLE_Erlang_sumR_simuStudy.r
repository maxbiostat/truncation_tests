library(sumR)
library(bbmle)
source("Erlang_aux.r")
#############################
nll_adaptive <- function(log_mu, log_beta, data){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood(data,
                          pars = theta,
                          eps =.Machine$double.eps,
                          verbose = FALSE)
}
nll_fixed <- function(log_mu, log_beta, data, Nit = 1000){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood(data,
                          adaptive = FALSE,
                          pars = theta,
                          eps =.Machine$double.eps,
                          N_fix = Nit,
                          verbose = FALSE)
}
nll_adaptive_full <- function(log_mu, log_beta, data){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full_fast(data,
                                    pars = theta,
                                    eps =.Machine$double.eps,
                                    verbose = FALSE)
}
nll_fixed_full <- function(log_mu, log_beta, data, Nit = 1000){
  theta <- exp(c(log_mu, log_beta))
  -marginal_loglikelihood_full_fast(data,
                                    adaptive = FALSE,
                                    pars = theta,
                                    eps =.Machine$double.eps,
                                    N_fix = Nit,
                                    verbose = FALSE)
}

get_estimates <- function(the.data, minusloglik_fun,
                          true.Mu = NA, true.Beta = NA){
  
  mod_nll <- function(log_mu, log_beta){
    minusloglik_fun(log_mu = log_mu,
                    log_beta = log_beta,
                    data = the.data)
  }
  
  inits <- log(
    c(runif(1, 1, 100),
      runif(1, .1, 50))
  )
  
  running.time <- system.time(
    MLE <- bbmle::mle2(minuslogl = mod_nll,
                       start = list(log_mu = inits[1],
                                    log_beta = inits[2]),
                       method = "L-BFGS-B"
    )
  )
  
  summy <- summary(MLE)
  
  point.ests <- as.numeric(summy@coef[, 1])
  MLE.se <- as.numeric(summy@coef[, 2])
  conf.level <- 0.95
  Z <- qnorm(p = (1 + conf.level)/2)
  
  MLE.lower <- point.ests - Z*MLE.se
  MLE.upper <- point.ests + Z*MLE.se
  
  MLE.ests <- tibble::tibble(
    parameter = c("mu", "beta"),
    true = c(true.Mu, true.Beta),
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
simulate_and_fit_once <- function(i,
                                  Mu,
                                  Beta,
                                  K){
  ## Keeping a log to spot problems
  logname <- paste0("Mu=", Mu, "_Beta=", Beta, "_K=", K, ".log")
  write(paste0("Doing replicate ", i),
        file = logname, append = TRUE)
  ##
  
  dt <- simulate_obsdata(n = K, mu = Mu, b = Beta,
                         seed = i)
  
  # fit.bessel.fixed <- get_estimates(the.data = dt,
  #                                   minusloglik_fun = nll_fixed,
  #                                   true.Mu = Mu, true.Beta = Beta)
  # 
  # fit.bessel.fixed$results$implementation <- "Bessel"
  # fit.bessel.fixed$results$approach <- "fixed"
  
  fit.bessel.adaptive <- get_estimates(the.data = dt,
                                       minusloglik_fun = nll_adaptive,
                                       true.Mu = Mu, true.Beta = Beta)

  fit.bessel.adaptive$results$implementation <- "Bessel"
  fit.bessel.adaptive$results$approach <- "adaptive"
  
  # fit.full.fixed <- get_estimates(the.data = dt,
  #                                 minusloglik_fun = nll_fixed_full,
  #                                 true.Mu = Mu, true.Beta = Beta)
  # 
  # fit.full.fixed$results$implementation <- "Full"
  # fit.full.fixed$results$approach <- "fixed"
  
  fit.full.adaptive <- get_estimates(the.data = dt,
                                     minusloglik_fun = nll_adaptive_full,
                                     true.Mu = Mu, true.Beta = Beta)

  fit.full.adaptive$results$implementation <- "Full"
  fit.full.adaptive$results$approach <- "adaptive"

  ## Assembly
  
  out <- do.call(rbind, 
                 list(
                   # fit.bessel.fixed$results,
                   fit.bessel.adaptive$results,
                   # fit.full.fixed$results,
                   fit.full.adaptive$results
                 ))
  out$replicate <- i
  return(out)
}
simulate_and_fit_reps <- function(Nrep, Mu, Beta, K, ncores = 2){
  
  logname <- paste0("Mu=", Mu, "_Beta=", Beta, "_K=", K, ".log")
  system(
    paste("rm ", logname) ## deleting old log files
  )
  
  raw <- parallel::mclapply(1:Nrep,
                            function(j){
                              simulate_and_fit_once(i = j,
                                                    Mu = Mu,
                                                    Beta = Beta,
                                                    K = K)
                            },
                            mc.cores = ncores)
  return(
    do.call(rbind, raw)
  )
}

Nreps <- 200
simus <- simulate_and_fit_reps(Nrep = Nreps,
                               Mu = 1500,
                               Beta = 1,
                               K = 5,
                               ncores = 8)
simus

# problematic <- setdiff(1:Nreps, unique(simus$replicate))
# 
# simulate_and_fit_once(
#   i = problematic[1],
#   Mu = 1500,
#   Beta = 1,
#   K = 5
# )
