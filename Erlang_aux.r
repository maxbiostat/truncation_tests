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
                     adaptive = TRUE, N_fix = 1000, verbose = FALSE){
  
  mu <- pars[1]
  b <- pars[2]
  log_z <- log(mu) + log(b) + log(x)

  if(adaptive){
    bessel.sumR <- sumR::infiniteSum(logFunction = "bessel_I_logX",
                                           parameters = c(log(2) + log_z/2, 1),
                                           epsilon = eps)
    if(verbose) cat("Marginalised probability took", bessel.sumR$n,
                    "iterations \n")
  }else{
    bessel.sumR <- sumR::finiteSum(logFunction = "bessel_I_logX",
                                         parameters = c(log(2) + log_z/2, 1),
                                         n = N_fix)
  }
  ans <- -(mu + b*x) - log(x) + log_z/2 + bessel.sumR$sum
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
                     adaptive = TRUE, N_fix = 100, verbose = FALSE){
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
    ans <- logProb$sum
  }
  return(ans - log1p(-exp(-pars[1])))
}

marginal_loglikelihood_full <- function(data, pars, 
                                   eps = .Machine$double.eps,
                                   adaptive = TRUE,
                                   N_fix = 1000,
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

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double erlang(long n, double *p)
{
  return n < 1 ? -INFINITY :
    R::dpois(n, p[0], 1) + R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double sum_erlang(Rcpp::NumericVector parameters, double epsilon = 1E-16,
  long maxIter = 100000)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = infiniteSum(erlang, parameter, -INFINITY, 0,
                  epsilon, maxIter, 0, &n);
  // Rcpp::Rcout << "Summation took " << n << " iterations to converge.\\n";

  return (double)r;
}
')

Rcpp::sourceCpp(code='
#include <Rcpp.h>

// [[Rcpp::depends(sumR)]]

#include <sumR.h>

long double erlang(long n, double *p)
{
  return n < 1 ? -INFINITY :
    R::dpois(n, p[0], 1) + R::dgamma(p[2], n, 1 / p[1], 1);
}

// [[Rcpp::export]]
double sum_erlang_fixed(Rcpp::NumericVector parameters, long iterations)
{
  double parameter[3];
  long double r;
  long n; // Number of iterations. Does not require initialization.

  parameter[0] = parameters[0];
  parameter[1] = parameters[1];
  parameter[2] = parameters[2];

  r = sumNTimes(erlang, parameter, iterations, 0);
  
  return (double)r;
}
')

marg_lik_full_fast <- function(x, pars, eps = .Machine$double.eps,
                          adaptive = TRUE, N_fix = 1000,
                          verbose = FALSE){
  if(adaptive){
    logProb <- sum_erlang(parameters = c(pars, x), epsilon = eps)
  }else{
    logProb <- sum_erlang_fixed(parameters = c(pars, x), N_fix)
  }
  return(logProb - log1p(-exp(-pars[1])))
}

marginal_loglikelihood_full_fast <- function(data, pars, 
                                        eps = .Machine$double.eps,
                                        adaptive = TRUE,
                                        N_fix = 1000,
                                        verbose = FALSE){
  logliks <- unlist(
    lapply(1:data$K, function(i){
      marg_lik_full_fast(x = data$obs_x[i], pars = pars,
                    eps = eps, adaptive = adaptive,
                    N_fix = N_fix, verbose = verbose) 
    })
  )
  
  return(sum(logliks))
}
