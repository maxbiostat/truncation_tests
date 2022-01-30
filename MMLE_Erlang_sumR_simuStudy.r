library(sumR)
library(bbmle)
source("Erlang_aux.r")
source("Erlang_hessian.R")
#############################
get_estimates <- function(data,
                          Mu, Beta,
                          fn = "BesselAdaptive",
                          conf.level = 0.95,
                          compute_Hessian = TRUE){
  
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
                            eps = .Machine$double.eps,
                            N_fix = 1000,
                            verbose = FALSE)
  }
  
  nll_fixed_full <- function(log_mu, log_beta){
    theta <- exp(c(log_mu, log_beta))
    -marginal_loglikelihood_full_fast(data,
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
  
  minusloglik_fun <- switch(fn,
    "BesselFixed" = nll_fixed,
    "BesselAdaptive" = nll_adaptive,
    "FullFixed" = nll_fixed_full,
    "FullAdaptive" = nll_adaptive_full_fast
  )
  
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
    true = c(Mu, Beta),
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
get_all_estimates <- function(dt, Mu, Beta){

  est.fixed <- get_estimates(data = dt,
                             Mu = Mu, Beta = Beta,
                             fn = "BesselFixed")
  est.adaptive <- get_estimates(data = dt,
                                Mu = Mu, Beta = Beta,
                                fn = "BesselAdaptive")
  est.fixed.full <- get_estimates(data = dt,
                                  Mu = Mu, Beta = Beta,
                                  fn = "FullFixed")
  est.adaptive.full.fast <- get_estimates(data = dt,
                                          Mu = Mu, Beta = Beta,
                                          fn = "FullAdaptive")
  outl <- list(
    tibble::tibble(est.fixed$results,
                   implementation = "Bessel",
                   method = "Fixed"),
    tibble::tibble(est.adaptive$results,
                   implementation = "Bessel",
                   method = "Adaptive"),
    tibble::tibble(est.fixed.full$results,
                   implementation = "Full",
                   method = "Fixed"),
    tibble::tibble(est.adaptive.full.fast$results,
                   implementation = "Full",
                   method = "Adaptive")
  )
  return(
    do.call(rbind, outl)
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
  
  dt <- simulate_obsdata(n = K, mu = Mu, b = Beta,
                         seed = i)
  out <- get_all_estimates(dt = dt,
                           Mu = Mu,
                           Beta = Beta)
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
  save(raw, 
       file = paste0("Mu=", Mu, "_Beta=", Beta, "_K=", K, ".RData"))
  return(
    do.call(rbind, raw)
  )
}

MM <- 1500
BB <- .1
J <- 50 
  
Nreps <- 500
simus <- simulate_and_fit_reps(Nrep = Nreps,
                               Mu = MM,
                               Beta = BB,
                               K = J,
                               ncores = 6)
write.csv(simus,
          file = paste0("Mu=", MM, "_Beta=", BB, "_K=", J, ".csv"),
          row.names = FALSE)