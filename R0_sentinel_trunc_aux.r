get_marginalised_prob <- function(r, w, p, verbose = FALSE){
  ## Get Pr(S = 0 | r, w, p)
  a <- r/w
  b <- r/(w * (1 + r/w)^(w + 1)) * (1-p)
  solve_for_u <- function(b, w){
    ## We know b < 1 because a and omega are positive
    ## it suffices to find the positive root smaller than 1/omega, which
    ## is the positive root of the objective function obj_fun
    b_of_u <- function(u) u/(1 + u)^{w + 1}
    obj_fun <- function(cand) (b-b_of_u(cand))^2
    u <- optimise(obj_fun, c(0, 1/w), tol = 1e-25)$minimum
    return(u)
  }
  ustar <- solve_for_u(b, w = w)
  if(verbose) cat("u* is:", ustar, "\n")
  return( (a + 1)/a * ustar/(1 + ustar))
}
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
log_diff_exp <- function(x, y) {
  # if(x <= y) stop("computing the log of a negative number")
  if(y == -Inf){
    return (x)
  }else{
    return (x + log1p(-exp(y-x)) )
  }
}
log1m_exp <- function(y){
  return(
    log_diff_exp(0, y)
  )
} 
obs.p.sentinel <- function(x, nu, log = FALSE){
  ans <- log1m_exp( x*log_diff_exp(0, log(nu)) )
  if(!log) ans <- exp(ans)
  return(ans)
} 
simulate_sentinel_data <- function(n, r, w, nu, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  Y <- rR0(n = n, r = r, w = w, D = 1)
  Z <- sapply(Y, function(x) rbinom(1, size = 1,
                                     prob = obs.p.sentinel(x = x, nu = nu)))
  X <- ifelse(Z == 1, Y, 0)
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

## Pr(Y = k | R0, omega) * Pr(X = x | Y = y, nu)
sentinel_summand <- function(k, theta){
  ans <- dR0(y = k, r = theta[1], w = theta[2], log = TRUE) + k*log1p(-theta[3])
  return(ans)
}
sentinel_lpmf <- function(k, theta){
  ans <- dR0(y = k, r = theta[1], w = theta[2], log = TRUE) + 
    log1m_exp(k * log1m_exp(log(theta[3])))
  return(ans)
}

## Pr(X = x | R0, omega, nu)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  
  logLL <- log(pars[1]) +
    (1 + pars[2]) * (log1p(pars[2]) - log(pars[1] + pars[2]) ) + log1p(-pars[3])
  
  if(adaptive){
    logPr0 <- sumR::infiniteSum(logFunction = sentinel_summand,
                                           parameters = c(pars, x),
                                           logL = logLL,
                                           epsilon = eps,
                                           n0 = 1,
                                           maxIter = 1E5)
    if(verbose) cat("Marginalised probability took", logPr0$n,
                    "iterations \n")
    
  }else{
    logPr0 <- sumR::finiteSum(logFunction = sentinel_summand,
                                         parameters = c(pars, x),
                                         n0 = 1,
                                         n = N_fix)
  }
  ans <- logPr0$sum
  return(ans)
}

marginal_loglikelihood_adaptive <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE, N_fix = NULL,
                                   verbose = FALSE){
  lp0 <- marg_lik(x = 0, pars = pars,
                 eps = eps, adaptive = adaptive,
                 N_fix = N_fix, verbose = verbose)
  
  logliks <- unlist(
    lapply(1:data$K, function(i){
      data$freq[i] * (sentinel_lpmf(k = i, theta = pars) - log1p(-exp(lp0)))
    })
  )
  return(sum(logliks))
}

marginal_loglikelihood_exact <- function(data, pars){
  p0 <- get_marginalised_prob(r = pars[1], w = pars[2], p = pars[3])
  lpp0 <- log1p(-p0)
  logliks <- unlist(
    lapply(1:data$K, function(i){
      data$freq[i] * ( sentinel_lpmf(k = i, theta = pars) - lpp0)
    })
  )
  logLikelihood <- sum(logliks)
  # cat("r:", pars[1], " w:", pars[2], " nu:", pars[3], "ll:", logLikelihood, "\n")
  return(logLikelihood)
}