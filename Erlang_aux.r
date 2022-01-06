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


Rcpp::cppFunction(code='
  NumericVector marginal_Erlang_lpdf_C(IntegerVector n, NumericVector p)
  {
  NumericVector output(n.size());
  for (int i = 0; i < n.size(); i++) {
    if (n[i] < 1) output[i] = -INFINITY;
    else
      output[i] = R::dpois(n[i], p[0], true) + R::dgamma(p[2], n[i], 1/p[1], true);
  }
  return output;
  }')

marg_lik_full <- function(x, pars, eps = .Machine$double.eps,
                     adaptive = TRUE, N_fix = NULL, verbose = FALSE){
  if(adaptive){
    logProb <- sumR::infiniteSum(logFunction = marginal_Erlang_lpdf_C,
                                     parameters = c(pars, x),
                                     logL = log(0),
                                     epsilon = eps)
    if(verbose) cat("Marginalised probability took", logProb$n,
                    "iterations \n")
    ans <- logProb$sum
  }else{
    logProb <- sumR::finiteSum(logFunction = marginal_Erlang_lpdf_C,
                                   parameters = c(pars, x),
                                   n = N_fix)
    ans <- logProb 
  }
  return(ans - log1p(-exp(-pars[1])))
}

marginal_loglikelihood_full <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE,
                                   N_fix = NULL,
                                   verbose = FALSE){
  logliks <- unlist(
    lapply(1:data$K, function(i){
      marg_lik_full(x = data$obs_x[i], pars = pars,
               eps = eps, adaptive = adaptive,
               N_fix = N_fix, verbose = verbose) 
    })
  )
  
  return(sum(logliks))
}