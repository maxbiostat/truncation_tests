library(sumR)
library(extraDistr)
library(rstan)
source("Poisson_Betabinomial_trunc_aux.r")
adaptive_model <- rstan::stan_model(file = "stan/Poisson_BetaBinomial_MLE_adaptive.stan",
                           allow_undefined = TRUE)
fixed_model <- rstan::stan_model(file = "stan/Poisson_BetaBinomial_MLE_fixed.stan",
                              allow_undefined = TRUE)

# cmdstanr::cmdstan_model("stan/Poisson_BetaBinomial_MLE_fixed.stan",
#                         include_paths = "./stan/")
get_p <- function(pars){
  pars[2]/(pars[2] + pars[3])
}

is_in <- function(x, l, u){
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}
####################
obtain_fit <- function(data,
                       theta,
                       target.eps,
                       M,
                       N_fix,
                       adaptive = TRUE){
    
  PBB.data <- list(
    K = data$K,
    x = data$count,
    n = data$freq,
    epsilon = target.eps,
    max_iter = M,
    fixed_nit = N_fix
  )
  
  if(adaptive){
    runtime <- system.time(
      Opt.MLE <- rstan::optimizing(obj = adaptive_model,
                                   data = PBB.data,
                                   hessian = FALSE,
                                   draws = 500,
                                   verbose = FALSE)
    )
    
  }else{
    runtime <- system.time(
      Opt.MLE <- rstan::optimizing(obj = fixed_model,
                                   data = PBB.data,
                                   hessian = FALSE,
                                   draws = 500,
                                   verbose = FALSE)
    )
  }
  
  point.ests <- as.numeric(Opt.MLE$par[1:3])
  conf.level <- 0.95
  
  draws.mean <- colMeans(Opt.MLE$theta_tilde)[1:3]
  draws.lwr <-  apply(Opt.MLE$theta_tilde, 2,
                      quantile, probs = (1 - conf.level)/2)[1:3]
  draws.upr <- apply(Opt.MLE$theta_tilde, 2,
                     quantile, probs = (1 + conf.level)/2)[1:3]
  
  draws.p <- apply(Opt.MLE$theta_tilde[, 1:3], 1,  get_p)
  
  pl <- quantile(draws.p, probs = (1 - conf.level)/2)
  pu <- quantile(draws.p, probs = (1 + conf.level)/2)  
  
  iters <- Opt.MLE$par[grep("niter", names(Opt.MLE$par))]
  
  MLE.ests <- tibble::tibble(
    parameter = c("mu", "alpha", "beta", "p_hat",
                  "log_lik", "n_iter"),
    true = c(theta[1], theta[2],
             theta[3], theta[2]/(theta[2]+theta[3]), NA,  NA),
    point =  c(point.ests, get_p(point.ests),
               Opt.MLE$value, median(iters)),
    lwr = c(draws.lwr, pl, NA, min(iters)),
    upr = c(draws.upr, pu, NA, max(iters)),
  )
  MLE.ests$covers <- is_in(x = MLE.ests$true, 
                           l = MLE.ests$lwr,
                           u = MLE.ests$upr)
  MLE.ests$time <- runtime[3]
  MLE.ests$fail <- as.logical(Opt.MLE$return_code)
  MLE.ests$method <- ifelse(adaptive, "adaptive", "fixed")
  return(MLE.ests)
}

simulate_and_fit_once <- function(true.mu,
                                  true.alpha,
                                  true.beta,
                                  N,
                                  target.eps = .Machine$double.eps,
                                  M = 2E4,
                                  N_fix = 1E3){
  simu <- simulate_obsdata(n = N,
                           mu = true.mu,
                           a = true.alpha,
                           b = true.beta)
  the.data <- compress_counts(simu$obs_xprime)
  
  fit.adapt <- obtain_fit(data = the.data,
                          theta = c(true.mu, true.alpha, true.beta),
                          target.eps = target.eps,
                          M = M,
                          N_fix = N_fix,
                          adaptive = TRUE)
  fit.fixed <- obtain_fit(data = the.data,
                          theta = c(true.mu, true.alpha, true.beta),
                          target.eps = target.eps,
                          M = M,
                          N_fix = N_fix,
                          adaptive = FALSE)
  
  return(rbind(fit.adapt, fit.fixed))
}
######

Nrep <- 200
results <- do.call(rbind, 
                   lapply(1:Nrep, function(i){
                     raw <- simulate_and_fit_once(
                       true.mu = 100,
                       true.alpha = 2,
                       true.beta = 2,
                       N = 500,
                       M = 2E4,
                       target.eps = .Machine$double.eps,
                       N_fix = 1E3
                     )
                     raw$replicate <- i
                     return(raw)
                   }))
results
write.csv(results, "../PBB_MMLE_C.csv", row.names = FALSE)
