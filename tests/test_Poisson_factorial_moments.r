library(sumR)
source("../aux/aux.r")
#######
poisson_lfactmom <- function(k, theta){
  x <- theta[2]
  ans <- ifelse(k < x, -Inf,
                lfactorial(k) - lfactorial(k-x) +
                  dpois(x = k, lambda = theta[1], log = TRUE))
    return(ans)
}
#######
Mu <- 12
R <- 2
Theta <- c(Mu, R)
TrueValue <- R*log(Mu)
Eps <- .Machine$double.eps
lgL <- log(0)
#######

result <- compare_approximations(compute_lterm = poisson_lfactmom,
                                 theta = Theta,
                                 exact = TrueValue,
                                 eps = Eps,
                                 logL = lgL)
result

