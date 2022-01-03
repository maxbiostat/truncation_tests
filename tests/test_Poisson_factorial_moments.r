library(sumR)
source("../aux/aux.r")
#######
poisson_lfactmom_R <- function(k, theta){
  x <- theta[2]
  ans <- ifelse(k < x, -Inf,
                lfactorial(k) - lfactorial(k-x) +
                  dpois(x = k, lambda = theta[1], log = TRUE))
  return(ans)
}

Rcpp::cppFunction(code='
  NumericVector poisson_lfactmom_C(IntegerVector n, NumericVector p)
  {
  NumericVector output(n.size());
  for (int i = 0; i < n.size(); i++) {
    if (p[1] > n[i]) output[i] = -INFINITY;
    else
      output[i] = R::dpois(n[i], p[0], true) +
      lgamma(n[i] + 1) - lgamma(n[i] - p[1] + 1);
  }
  return output;
  }')

#######
Mu <- 11.8 ## try 11.8 and 11.9 and 12.8 and be amazed
R <- 2
Theta <- c(Mu, R)
TrueValue <- R*log(Mu)
Eps <- .Machine$double.eps*1E4
lgL <- log(0)
#######

result.R <- compare_approximations(compute_lterm = poisson_lfactmom_R,
                                 theta = Theta,
                                 exact = TrueValue,
                                 eps = Eps,
                                 logL = lgL)

result.C <- compare_approximations(compute_lterm = poisson_lfactmom_C,
                                   theta = Theta,
                                   exact = TrueValue,
                                   eps = Eps,
                                   logL = lgL)


result.preComp <- compare_approximations(
                       compute_lterm = "poisson_fact_moment",
                       theta = Theta,
                       exact = TrueValue,
                       eps = Eps)

result.R
result.C
result.preComp

####
##$ Table
spit_approx <- function(lambda, order, epsilon){
  ans <- suppressWarnings(
    compare_approximations(
      compute_lterm = "poisson_fact_moment",
      theta = c(lambda, order),
      exact = order*log(lambda),
      max_iter = 3E5,
      eps = epsilon)
  ) 
  out <- tibble::tibble(
    mu = lambda,
    r = order,
    target_error = ans$target_error,
    error = ans$error,
    success = abs(ans$error) <= ans$target_error,
    n_iter = ans$n_evaluations,
    method = ans$Method
  )
  return(out)
}

approximation.grid <- expand.grid(
  Mu = c(0.5, 1, 10, 100),
  R = c(2, 3, 5, 10),
  eps = 10^(0:6)*.Machine$double.eps
)

approximation.list <- lapply(
  1:nrow(approximation.grid),
  function(i){
    linha <- approximation.grid[i, ] ## linha == Portuguese for row
    spit_approx(lambda = linha$Mu,
                order = linha$R,
                epsilon = linha$eps)      
  })

approximation.dt <- do.call(rbind, approximation.list)

aggregate(success~method, approximation.dt, mean)
aggregate(success~method+target_error, approximation.dt, mean)


#### Investigating the instances in which adaptive/naive fails
#### but the sum with a huge (3E5) number of terms does not.

approximation.dt$row_id <- paste(approximation.dt$mu,
                                  approximation.dt$r,
                                  approximation.dt$target_error,
                                  sep = "_")
( fail.adapt <- subset(approximation.dt, method == "Adaptive" & !success) )
( fail.naive <- subset(approximation.dt, method == "Naive" & !success) )
( fail.huge <- subset(approximation.dt, method == "Fixed_3e+05" & !success) )

interesting.naive <- subset(approximation.dt,
                            row_id %in% setdiff(fail.naive$row_id, fail.huge$row_id))

interesting.adaptive <- subset(approximation.dt,
                               row_id %in% setdiff(fail.adapt$row_id, fail.huge$row_id))

interesting.adaptive
tail(interesting.adaptive)
