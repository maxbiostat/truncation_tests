PL_diff_function_direct <- function(k, theta){
  alpha <- theta[1]
  xmin <- theta[2]
  phi <- theta[3:4]
  ans <- sapply(k, function(k) ifelse(k < xmin,  -Inf,
                                      -alpha*log(k) + logDiffExp(0, -phi[1] - phi[2] * (k-xmin)))
  )
  return(ans)
}