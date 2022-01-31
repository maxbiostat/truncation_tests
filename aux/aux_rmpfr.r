log1m <- function(x) {
  Rmpfr::log1mexp(log(x))
}
##
log1m_exp <- function(a) {
  Rmpfr::log1mexp(a)
}
##
log_diff_exp <- function(x, y){
  M <- max(x, y)
  m <- min(x, y)
  A <- M
  B <- m-M
  return(
    A + Rmpfr::log1mexp(-B)
  )
}
##
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- log_diff_exp(x, y)
  }else{
    ans <- sgn * exp(log_diff_exp(x, y))
  }
  return(ans)
}
##
relative_difference <- function(x, y){
  obj <- Rmpfr::all.equal(x, y, tolerance = 0)
  ans <- suppressWarnings(
    Rmpfr::mpfr(as.numeric(gsub("Mean relative difference: ", "", obj)), 1024)
  )
  if(is.na(ans)) ans <- Rmpfr::mpfr(0, 1024)
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
  Sums <- c(Rmpfr::mpfr(naive$sum, 1024),
            # batches$sum,
            Rmpfr::mpfr(batches_C$sum, 1024),
            Rmpfr::mpfr(adaptive$sum, 1024),
            Rmpfr::mpfr(fixed$sum, 1024),
            Rmpfr::mpfr(fixedMax$sum, 1024))
  Niters <- c(naive$n,
              # batches$n,
              batches_C$n,
              adaptive$n,
              Nfix,
              max_iter)
  K <- length(Sums)
  ##
  out <- list(
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
    target_error = rep(mpfr(eps, 1024), K),
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
