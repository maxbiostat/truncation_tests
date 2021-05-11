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
                                   epsilon = 1E-10,
                                   n0 = 0, max_iter = 1E5, logL = -Inf){
  naive <- adaptiveSum::naive_sum(lFun = compute_lterm,
                                  params = theta,
                                  eps = epsilon, maxIter = max_iter,
                                  n0 = n0) 
  adaptive <- adaptiveSum::adapt_sum(lFun = compute_lterm,
                                     params = theta,
                                     eps = epsilon, maxIter = max_iter,
                                     logL = logL, n0 = n0) 
  ##
  Sums <- c(naive$sum, adaptive$sum)
  Niters <- c(naive$n, adaptive$iter)
  ##
  out <- data.frame(
    Method = c("Naive", "Adaptive"),
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