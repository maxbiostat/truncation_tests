##
log1m <- function(x) {
  return(log1p(-x))
}
##
log1m_exp <- function(a) {
  if (a > -0.693147) {
    return(log(-expm1(a))) 
  } else {
    return(log1m(exp(a)));
  }
}
##
log_diff_exp <- function (x, y){
  return(x + log1m_exp(y - x))
} 
##
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- log_diff_exp(max(x, y), min(x, y))
  }else{
    ans <- sgn * exp(log_diff_exp(max(x, y), min(x, y)))
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
##
compare_approximations <- function(compute_lterm, theta, exact,
                                   eps = 1E-10,
                                   n0 = 0, Nstart = 10, Nfix = 1000,
                                   max_iter = 1E5, logL = -Inf){
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
  
  doubling <- sumR::infiniteSum_cFolding(
    logFunction = compute_lterm,
    parameters = theta,
    epsilon = eps,
    N_start = Nstart,
    c = 2,
    maxIter = max_iter,
    n0 = n0
  )
  
  doubling_C <- sumR::infiniteSum_cFolding_C(
    logFunction = compute_lterm,
    parameters = theta,
    epsilon = eps,
    N_start = Nstart,
    c = 2,
    maxIter = max_iter,
    n0 = n0
  )
  
  fixed <- sumR::finiteSum(
    logFunction = compute_lterm,
    parameters = theta,
    n0 = n0,
    n = Nfix
  )

  ##
  Sums <- c(naive$sum, doubling$sum, doubling_C$sum, adaptive$sum, fixed)
  Niters <- c(naive$n, doubling$n, doubling_C$n, adaptive$n, Nfix)
  ##
  out <- data.frame(
    Method = c("Naive", "C-folding", "C-folding_C", "Adaptive",
               paste0("Fixed_", Nfix)),
    error_logspace = Sums-exact,
    # relative_error_logspace = sapply(Sums,
    #                                  function(x) relative_difference(x, exact)),
    # relative_error = sapply(exp(Sums),
    #                         function(x) relative_difference(x, exp(exact))),
    log_error = sapply(Sums,
                       function(x) robust_difference(x, exact, log = TRUE)),
    error = sapply(Sums,
                   function(x) robust_difference(x, exact)),
    n_evaluations = Niters,
    sum = Sums
  )
  return(out)
}