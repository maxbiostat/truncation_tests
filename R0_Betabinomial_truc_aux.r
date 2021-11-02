dR0 <- function(y, r, w, log = FALSE){
  c1 = lgamma(w*y + y-1)
  c2 = lgamma(w*y)
  c3 = lgamma( y+1)
  c4 = (y-1) * (log(r) - log(w))
  c5 = (w*y + y -1) * log(1 + r/w)
  dens = c1 - (c2+c3) + (c4-c5)
  if(!log){
    dens <- exp(dens)
  }
  return(dens)
}
#
rR0 <- function(n, r, w, D, UpperBound = 1e4){
  ys <- 1:UpperBound
  Ps <- dR0(y = ys, r = r, w = w)
  as.vector(
    matrix(
      sample(ys, n * D, prob = Ps, replace = TRUE),
      ncol = D, nrow = n
    )  
  )
}
simulate_obsdata <- function(n, r, w, a, b, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  Y <- rR0(n = n, r = r, w = w, D = 1)
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
  ans <- dR0(y = k, r = theta[1], w = theta[2], log = TRUE) +
    extraDistr::dbbinom(x = theta[5], size = k,
                        alpha = theta[3], beta = theta[4], log = TRUE)
  
  return(ans)
}

## Pr(X = x | mu, a, p)
### pars = c(mu, a, p)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  
  logLL <- log(pars[1]) +
    (1 + pars[2]) * (log1p(pars[2]) - log(pars[1] + pars[2]) )
  
  if(adaptive){
    logMargProbKernel <- sumR::infiniteSum(logFunction = with_obs_error_lpmf,
                                           parameters = c(pars, x),
                                           logL = logLL,
                                           epsilon = eps,
                                           n0 = 1,
                                           maxIter = 1E5)
    if(verbose) cat("Marginalised probability took", logMargProbKernel$n,
                    "iterations \n")
    
  }else{
    logMargProbKernel <- sumR::finiteSum(logFunction = with_obs_error_lpmf,
                                         parameters = c(pars, x),
                                         n0 = 1,
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