source("truncation_package.r")
source("truncation_aux.r")

COMP_lpdf <- function(k, theta){
  lambda <- theta[1]
  nu <- theta[2]
  return(
    k * log(lambda) - nu*lfactorial(k)  
  )
}

Lambda <- .5
Nu <- .2
Theta <- c(Lambda, Nu)
Eps <- 1E-15
M <- 2E5
if(Nu == 1){
  TrueValue <- Lambda  
}else{
  if(Nu == 2){
    TrueValue <- log(besselI(2*sqrt(Lambda), nu = 0))
  }else{
    lps <- COMP_lpdf(k = 0:M, theta = Theta)
    TrueValue <- log_sum_exp(lps)
  }
}

TrueValue
result <- compare_approximations(COMP_lpdf, theta = Theta,
                                 exact = TrueValue, epsilon = Eps, max_iter = M/2)

result

# library(microbenchmark)
# 
# timing <- microbenchmark(
#   naive = approx_naive(COMP_lpdf, theta = Theta, epsilon = Eps, max_iter = M/2),
#   naive_threshold = approx_naive_under(COMP_lpdf, theta = Theta, epsilon = Eps, max_iter = M/2),
#   doubling = approx_doubling(COMP_lpdf, theta = Theta, epsilon = Eps, max_iter = M/2),
#   adaptive = approx_adaptive(COMP_lpdf, theta = Theta, epsilon = Eps, max_iter = M/2)
# )
# 
# result
# timing
