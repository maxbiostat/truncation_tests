simulate_obsdata <- function(n, mu, b, seed = NULL){
  if(!is.null(seed)) set.seed(seed)

  Y <- qpois(runif(n, dpois(0, mu), 1), mu)
  X <- rep(NA, n)
  
  for (i in 1:n){
    Zs <- rexp(n = Y[i], rate = b)
    X[i] <- sum(Zs)
  }
  
  return(
    list(
      obs_y = Y,
      obs_x = X,
      K = n
    )
  )
}


## Pr(X = x | mu, b, x)
### pars = c(mu, b, x)
marg_lik <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  
  mu <- pars[1]
  b <- pars[2]
  z <- mu*b*x

  
  if(adaptive){
    bessel.sumR <- sumR::infiniteSum(logFunction = "bessel_I",
                                           parameters = c(2*sqrt(z), 1),
                                           epsilon = eps)
    if(verbose) cat("Marginalised probability took", bessel.sumR$n,
                    "iterations \n")
    ans <- -(mu + b*x) - log(x) + log(z)/2 + bessel.sumR$sum  
  }else{
    bessel.sumR <- sumR::finiteSum(logFunction = "bessel_I",
                                         parameters = c(2*sqrt(z), 1),
                                         n = N_fix)
    ans <- -(mu + b*x) - log(x) + log(z)/2 + bessel.sumR
  }
  return(ans-log1p(-exp(-mu)))
}

marginal_loglikelihood <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE, N_fix = NULL,
                                   verbose = FALSE){
  logliks <- unlist(
    lapply(1:data$K, function(i){
      marg_lik(x = data$obs_x[i], pars = pars,
                              eps = eps, adaptive = adaptive,
                              N_fix = N_fix, verbose = verbose) 
    })
  )
  
  return(sum(logliks))
}