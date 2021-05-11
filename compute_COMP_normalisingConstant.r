library(adaptiveSum)
source("compare_approximations.r")
######################################
COMP_lpmf <- function(k, theta){
  lambda <- theta[1]
  nu <- theta[2]
  return(
    k * log(lambda) - nu*lfactorial(k)  
  )
}

Lambda <- 2
Nu <- .5
Theta <- c(Lambda, Nu)
Eps <- 1E-16
M <- 2E5
x0 <- 0
lgL <- log(0)

if(Nu == 1){
  TrueValue <- Lambda
}else{
  if(Nu == 2){
    TrueValue <- log(besselI(2*sqrt(Lambda), nu = 0))
  }else{
    lps <- COMP_lpmf(k = M:0, theta = Theta)
    TrueValue <- matrixStats::logSumExp(lps)
  }
}

TrueValue
result <- compare_approximations(COMP_lpmf, theta = Theta,
                                 exact = TrueValue,
                                 n0 = x0, logL = lgL,
                                 epsilon = Eps, max_iter = M/2)

result

adaptiveSum::naive_sum(lFun = COMP_lpmf,
                       params = Theta,
                       eps = Eps, maxIter = M/2,
                       n0 = x0)

timing <- bench::mark(
  naive = adaptiveSum::naive_sum(lFun = COMP_lpmf,
                                 params = Theta,
                                 eps = Eps, maxIter = M/2,
                                 n0 = x0),
  adaptive = adaptiveSum::adapt_sum(lFun = COMP_lpmf,
                                    params = Theta,
                                    eps = Eps, maxIter = M/2,
                                    n0 = x0, logL = lgL),
  check = FALSE
)

timing