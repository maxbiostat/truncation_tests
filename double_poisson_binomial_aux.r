simulate_obsdata <- function(n, mu, phi, p){
  require(gamlss.dist)
  Y <- rDPO(n = n, mu = mu, sigma = phi)
  X <- rbinom(n = n, size = Y, prob = p)
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


# theta = c(mu, phi, p)
# Pr(Y = k | mu, phi)
double_poisson_lpmf <- function(k, pars){
  mu <- pars[1]
  phi <- pars[2]
  lct <- log(phi)/2  - phi*mu
  ans <- sapply(k, function(k){
    ifelse(k == 0,  lct,
           lct -k + k*log(k) - lfactorial(k) + phi*k * (log(exp(1)*mu)-log(k)))
    
  }   )
  return(ans)
}

## Pr(Y = k | mu, phi) * Pr(X = x | Y = y, p)
### theta = c(mu, phi, p, x)
with_obs_error_lpmf <- function(k, theta){
  ans <- double_poisson_lpmf(k = k, pars = theta) +
    dbinom(x = theta[4], size = k, prob = theta[3], log = TRUE)
  return(ans)
}

## Pr(X = x | mu, phi, p)
### pars = c(mu, phi, p)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  if(adaptive){
    logConstDPO <- sumR::infiniteSum(logFunction = "double_poisson",
                                     parameters = pars[1:2],
                                     epsilon = eps,
                                     maxIter = 1E5)
    if(verbose) cat("Normalising constant took", logConstDPO$n,
                    "iterations \n")
    logMargProbKernel <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                                           parameters = c(pars, x),
                                           logL = log(0),
                                           epsilon = eps,
                                           maxIter = 1E5)
    if(verbose) cat("Marginalised probability took", logMargProbKernel$n,
                    "iterations \n")
    
  }else{
    logConstDPO <- sumR::finiteSum(logFunction = "double_poisson",
                                   parameters = pars[1:2],
                                   n = N_fix)
    logMargProbKernel <- sumR::finiteSum(logFunction = with_obs_error_lpmf,
                                         parameters = c(pars, x),
                                         n = N_fix)
  }
  ans <- logMargProbKernel$sum - logConstDPO$sum
  return(ans)
}

marginal_loglikelihood <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE, N_fix = NULL,
                                   verbose = FALSE){
  
  logliks <- unlist(
    lapply(1:data$K, function(i){
      data$freq[i] * marg_lik(x = data$count[i], pars = pars,
                              eps = eps, adaptive = adaptive,
                              N_fix = N_fix, verbose = verbose)
    })
  )
  return(sum(logliks))
}