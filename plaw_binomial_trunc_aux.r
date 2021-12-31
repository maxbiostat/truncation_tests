simulate_obsdata <- function(n, x_m, alpha, p, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  require(extraDistr)
  Y <- poweRlaw::rpldis(n = n, xmin = x_m, alpha = alpha)
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

plaw_ulpmf <- function(k, theta){
  if(theta[2] <= 1) return(rep(log(0), length(k)))
  lconst <- log(VGAM::zeta(theta[2], shift = theta[1]))
  lans <- ifelse(
    k > theta[1],
    yes = -theta[2]*log(k),
    no = log(0)
  ) - lconst
  return(lans)
}
## Pr(Y = k | x_m, alpha) * Pr(X = x | Y = y, p)
### theta = c(x_m, alpha, p, x)
with_obs_error_lpmf <- function(k, theta){
  ans <- plaw_ulpmf(k, theta) +
    dbinom(x = theta[4], size = k, prob = theta[3], log = TRUE)
  return(ans)
}

## Pr(X = x | x_m, alpha, p, x)
### pars = c(x_m, alpha, p, x)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  if(adaptive){
    logMargProbKernel <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                                           parameters = c(pars, x),
                                           logL = log1p(-pars[3]),
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
  
  # cat("parameters:", pars, " lp0:", lp0, "\n")
  
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