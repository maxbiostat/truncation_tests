source("truncation_package.r")
source("truncation_aux.r")

weird_series <- function(k, theta){
  eta <- theta[1]
  if(k==0) return(-Inf)
  return(
   lfactorial(k)-k*log(k)
  )
}
weird_series <- Vectorize(weird_series)

Theta <- NA
Eps <- 1E-20
TrueValue <- log_sum_exp(weird_series(0:1000, NA))#log(1.87985)

result <- compare_approximations(weird_series, theta = Theta,
                                 exact = TrueValue, epsilon = Eps, max_iter = 1E5)
result