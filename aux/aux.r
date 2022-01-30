Rcpp::cppFunction('double logDiffExp(double x, double y)
{return x > y ? 
Rf_logspace_sub(x, y) :
Rf_logspace_sub(y, x);}') 
##
log1m <- function(x) {
  return(log1p(-x))
}
##
log1m_exp <- function(a) {
  logDiffExp(0, a)
}
##
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- logDiffExp(x, y)
  }else{
    ans <- sgn * exp(logDiffExp(x, y))
  }
  return(ans)
}
##
relative_difference <- function(x, y){
  obj <- all.equal(x, y, tolerance = 0)
  ans <- suppressWarnings(
    as.numeric(gsub("Mean relative difference: ", "", obj))
  )
  if(is.na(ans)) ans <- 0
  return(ans)
}

compare_approximations <- function(compute_lterm, theta, exact,
                                   eps = 1E-10,
                                   n0 = 0,
                                   batch_size = 20,
                                   Nfix = 1000,
                                   max_iter = 3E5,
                                   logL = -Inf){
  naive <- sumR::infiniteSum(
    logFunction = compute_lterm,
    parameters = theta,
    epsilon = eps,
    maxIter = max_iter,
    n0 = n0,
    forceAlgorithm = 1
  )

  adaptive <- sumR::infiniteSum(
    logFunction = compute_lterm,
    parameters = theta,
    epsilon = eps,
    maxIter = max_iter,
    logL = logL,
    n0 = n0,
    forceAlgorithm = 2
  )

  # batches <- sumR::infiniteSum_batches(
  #   logFunction = compute_lterm,
  #   parameters = theta,
  #   epsilon = eps,
  #   batch_size = batch_size,
  #   maxIter = max_iter,
  #   n0 = n0
  # )

  batches_C <- sumR::infiniteSum_batches_C(
    logFunction = compute_lterm,
    parameters = theta,
    epsilon = eps,
    batch_size = batch_size,
    maxIter = max_iter,
    n0 = n0
  )

  fixed <- sumR::finiteSum(
    logFunction = compute_lterm,
    parameters = theta,
    n0 = n0,
    n = Nfix
  )

  fixedMax <- sumR::finiteSum(
    logFunction = compute_lterm,
    parameters = theta,
    n0 = n0,
    n = max_iter
  )

  ##
  Sums <- c(naive$sum,
            # batches$sum,
            batches_C$sum,
            adaptive$sum,
            fixed$sum,
            fixedMax$sum)
  Niters <- c(naive$n,
              # batches$n,
              batches_C$n,
              adaptive$n,
              Nfix,
              max_iter)
  K <- length(Sums)
  ##
  out <- data.frame(
    Method = c("Threshold",
               # "Batches",
               "Batches_C",
               "BoundingPair",
               paste0("Fixed_", Nfix),
               paste0("Fixed_", max_iter)),
    error_logspace = Sums-exact,
    # relative_error_logspace = sapply(Sums,
    #                                  function(x) relative_difference(x, exact)),
    percent_error = sapply(exp(Sums),
                           function(x) 100*relative_difference(x, exp(exact))),
    log_error = sapply(Sums,
                       function(x) robust_difference(x, exact, log = TRUE)),
    target_error = rep(eps, K),
    error = sapply(Sums,
                   function(x) robust_difference(x, exact)),
    error_max = sapply(Sums,
                   function(x) robust_difference(x, fixedMax$sum)),
    n_evaluations = Niters,
    sum = Sums,
    true = rep(exact, K)
  )
  return(out)
}
