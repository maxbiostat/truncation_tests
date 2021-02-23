##############
### Auxiliary functions
compare_approximations <- function(compute_lterm, theta, exact, epsilon = 1E-10,
                                   n0 = 0, max_iter = 1E5, logL = -Inf){
  naive <- approx_naive(compute_lterm, theta = theta, epsilon = epsilon, n0 = n0,
                        max_iter = max_iter)
  naive.u <- approx_naive_under(compute_lterm, theta = theta, epsilon = epsilon,
                                n0 = n0, max_iter = max_iter)
  doubling <- approx_doubling(compute_lterm, theta = theta, epsilon = epsilon,
                              n0 = n0, max_iter = max_iter)
  adaptive <- approx_adaptive(compute_lterm, theta = theta, epsilon = epsilon,
                              n0 = n0, max_iter = max_iter, logL = logL)
  ##
  Sums <- c(naive$sum, naive.u$sum, doubling$sum, adaptive$sum)
  Niters <- c(naive$iter, naive.u$iter, doubling$iter, adaptive$iter)
  ##
  out <- data.frame(
    Method = c("Naive", "Naive_threshold", "Doubling", "Adaptive"),
    absolute_error_logspace = Sums-exact,
    relative_error_logspace = sapply(Sums, function(x) relative_difference(x, exact)),
    relative_error = sapply(exp(Sums), function(x) relative_difference(x, exp(exact))),
    absolute_error = sapply(Sums, function(x) robust_difference(x, exact)),
    n_evaluations = Niters,
    sum = Sums
  )
  return(out)
}