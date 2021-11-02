simulate_obsdata <- function(n, mu, a, b, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  require(extraDistr)
  Y <- rpois(n = n, lambda = mu)
  X <- extraDistr::rbbinom(n = n, size = Y, alpha = a, beta = b)
  return(
    list(
      obs_y = Y,
      obs_x = X,
      obs_xprime = X[X>0]
    )
  )
}

compress_counts <- function(x){
  tt <- table(x)
  return(
    list(
      freq = as.numeric(tt),
      count = as.numeric(names(tt)),
      K = length(tt)
    )
  )
}

## Pr(Y = k | mu) * Pr(X = x | Y = y, a, b)
### theta = c(mu, a, b, x)
with_obs_error_lpmf <- function(k, theta){
  ans <-  dpois(x = k, lambda = theta[1], log = TRUE) +
    extraDistr::dbbinom(x = theta[4], size = k,
                        alpha = theta[2], beta = theta[3], log = TRUE)
  
  return(ans)
}

## Pr(X = x | mu, a, p)
### pars = c(mu, a, p)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  if(adaptive){
    logMargProbKernel <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                                           parameters = c(pars, x),
                                           logL = log(0),
                                           epsilon = eps,
                                           maxIter = 1E5)
    if(verbose) cat("Marginalised probability took", logMargProbKernel$n,
                    "iterations \n")
    
  }else{
    logMargProbKernel <- sumR::finiteSum(logFunction = with_obs_error_lpmf,
                                         parameters = c(pars, x),
                                         n = N_fix)
  }
  ans <- logMargProbKernel$sum
  return(ans)
}

marginal_loglikelihood <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE, N_fix = NULL,
                                   verbose = FALSE){
  lp0 <- marg_lik(x = 0, pars = pars,
                 eps = eps, adaptive = adaptive,
                 N_fix = N_fix, verbose = verbose)
  
  logliks <- unlist(
    lapply(1:data$K, function(i){
      data$freq[i] * (marg_lik(x = data$count[i], pars = pars,
                              eps = eps, adaptive = adaptive,
                              N_fix = N_fix, verbose = verbose) -
                        log1p(-exp(lp0)))
    })
  )
  return(sum(logliks))
}