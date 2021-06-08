# remotes::install_github("GuidoAMoreira/sumR")
library(sumR)
source("compare_approximations.r")
######################################
COMP_lpmf <- function(k, theta){
  lambda <- theta[1]
  nu <- theta[2]
  return(
    k * log(lambda) - nu*lfactorial(k)  
  )
}

Lambda <- 4
Nu <- .5
Theta <- c(Lambda, Nu)
Eps <- 1E-16
M <- 2E5
x0 <- 0
lgL <- -Inf 
Ns <- 10

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
result <- compare_approximations(compute_lterm = COMP_lpmf,
                                 theta = Theta,
                                 exact = TrueValue,
                                 n0 = x0, logL = lgL,
                                 Nstart = Ns,
                                 eps = Eps,
                                 max_iter = M/2)
result